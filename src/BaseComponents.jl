#
# Copyright (c) 2025 Andreas Hofmann
# Licensed under the MIT license. See LICENSE file in the project root for details.
#

# Basic components required for ScaraRobot example
# partially based on the Modelica Standard Library https://github.com/modelica/ModelicaStandardLibrary

module BaseComponents

    export EMF, ViscousFriction, SignalDCMotor


    using ModelingToolkit
    using ModelingToolkit: t_nounits as t, D_nounits as D
    using ModelingToolkitStandardLibrary
    using ModelingToolkitStandardLibrary.Electrical
    using ModelingToolkitStandardLibrary.Blocks
    using ModelingToolkitStandardLibrary.Mechanical: TranslationalModelica as Translational, Rotational



    @mtkmodel EMF begin
        @extend v,i = oneport = Electrical.OnePort()	
        @parameters begin
            k, [description = "Transformation coefficient"]
        end
        @variables begin
            phi(t), [description = "Angle"]
            w(t), [description = "Angular Velocity"]	
        end
        @components begin
            flange = Rotational.Flange()
        end	
        @equations begin
            phi ~ flange.phi
            D(phi) ~ w
            k * w ~ v
            flange.tau ~ -k * i	
        end
    end

    @mtkmodel ViscousFriction begin
        @parameters begin
            frictionFactor = 6.8e-5, [description = "Angular Velocity dependent Friction Factor [Nm.s/rad] "]
        end
        @components begin
            flange = Rotational.Flange()
        end
        @equations begin
            flange.tau ~ frictionFactor * D(flange.phi) # since friction force is negatively applied it acts inside the component rather than outside. Therefore, no negative sign is needed
        end
    end


    @mtkmodel SignalDCMotor begin
        @parameters begin
            R = 0.156, [description = "Motor Resistance [Ω]"]
            L = 95e-6, [description = "Motor Inductivity [H]"]
            k = 0.02, [description = "Motor Constant [N.m/A]"]
            J = 2.86e-5, [description = "Motor Inertia [kg.m2]"]
            v = 6.8e-5, [description = "Viscous friction factor [Nm.s/rad]"]
        end
        @components begin
            vSignal = Blocks.RealInput()
            signalVoltage = Electrical.Voltage()
            resistor = Electrical.Resistor(;R = R)
            inductor = Electrical.Inductor(;L=L)
            emf = BaseComponents.EMF(;k = k)
            ground = Electrical.Ground()  
            currentSensor = Electrical.CurrentSensor()
            inertia = Rotational.Inertia(;J=J)
            flange_a = Rotational.Flange()
            viscousFriction = BaseComponents.ViscousFriction(;frictionFactor=v)
        end
        @equations begin
            connect(vSignal, signalVoltage.V)
            connect(signalVoltage.p, resistor.p)
            connect(resistor.n, inductor.p)
            connect(inductor.n, emf.p)
            connect(emf.n, currentSensor.p)
            connect(currentSensor.n, signalVoltage.n, ground.g)
            connect(emf.flange,inertia.flange_a)
            connect(inertia.flange_b,flange_a)
            connect(inertia.flange_a,viscousFriction.flange)
        end
    end


end