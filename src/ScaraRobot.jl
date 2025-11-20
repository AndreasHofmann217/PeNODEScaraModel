#
# Copyright (c) 2025 Andreas Hofmann
# Licensed under the MIT license. See LICENSE file in the project root for details.
#


module ScaraRobot
    using ModelingToolkit
    using ModelingToolkit: t_nounits as t, D_nounits as D
    using ModelingToolkitStandardLibrary.Blocks
    using ModelingToolkitStandardLibrary.Mechanical: TranslationalModelica as Translational, Rotational

    using NaNMath

    
    using ..BaseComponents
    using ..PartialPlanarMechanics
    using ..SymbolicInputFunctions
    using ..ScaraKinematic


    export ScaraMTKModel, FrictionlessScaraMTKModel




    "Friction based TCP"
    @mtkmodel TCP begin
        @parameters begin
            eps = 1e-6
            mu_A = 0.3
            mu_S = 0.1
            vSlide = 0.001
            delta = 2
            viscousFrictionFactor = 5e-4
        end
        @variables begin
            x(t)
            y(t)
            vx(t)
            vy(t)
            f_abs(t)
            v_abs(t)
            fx(t)
            fy(t)
            Fs(t)
            Fc(t)
            Fv(t)
            Fstribeck(t)
        end
        @components begin
            frame_a = PartialPlanarMechanics.Frame()
            zForceIn = Blocks.RealInput()
            posXOut = Blocks.RealOutput()
            posYOut = Blocks.RealOutput()
        end
        @equations begin
            frame_a.tau ~ 0.0
            frame_a.fx ~ fx
            frame_a.fy ~ fy
            
            frame_a.x ~ posXOut.u
            frame_a.y ~ posYOut.u
    
            frame_a.x ~ x
            frame_a.y ~ y
            D(x) ~ vx
            D(y) ~ vy
    
            # implementation of stribeck friction based upon
            # https://mogi.bme.hu/TAMOP/robot_applications/ch07.html#ch-8.3.4
            Fc ~ zForceIn.u * mu_S
            Fs ~ zForceIn.u * mu_A
            Fv ~ viscousFrictionFactor * v_abs
            Fstribeck ~ Fc + (Fs-Fc)*exp(-(v_abs/vSlide)^delta) + Fv
            f_abs ~ Fstribeck 
            v_abs ~ NaNMath.sqrt(vx^2 + vy^2)
            fx ~ ifelse(v_abs > eps,-f_abs * vx/v_abs,0.0)
            fy ~ ifelse(v_abs > eps,-f_abs * vy/v_abs,0.0)
        end
    end


    "Friction less TCP"
    @mtkmodel FrictionlessTCP begin
        @parameters begin
        end
        @variables begin
            x(t)
            y(t)
            vx(t)
            vy(t)
            fx(t)
            fy(t)
        end
        @components begin
            frame_a = PartialPlanarMechanics.Frame()
            zForceIn = Blocks.RealInput()
            posXOut = Blocks.RealOutput()
            posYOut = Blocks.RealOutput()
        end
        @equations begin
            frame_a.tau ~ 0.0
            frame_a.fx ~ fx
            frame_a.fy ~ fy
            
            frame_a.x ~ posXOut.u
            frame_a.y ~ posYOut.u
    
            frame_a.x ~ x
            frame_a.y ~ y
            D(x) ~ vx
            D(y) ~ vy
            fx ~ 0.0
            fy ~ 0.0
        end
    end


    @mtkmodel ScaraMechanics begin
        @parameters begin
            m1 = 0.008, [description="mass of Link1"]
            J1 = 2.7e-5, [description="Inertia of Link1"]
            L1 = 0.2, [description="Length of Link1"]
            m2 = 0.004, [description="mass of Link2"]
            J2 = 3.5e-6, [description="Inertia of Link2"]
            L2 = 0.1, [description="Length of Link2"]
            eps = 1e-6, [description="Stiction threshold"]
            vSlide = 0.001, [description = "Sliding velocity"]
            mu_A = 0.3, [description = "Friction coefficient at adhesion"]
            mu_S = 0.1, [description = "Friction coefficient at sliding"]
        end
        @components begin
            fixed = PartialPlanarMechanics.Fixed()
            rev1 = PartialPlanarMechanics.RevoluteJointActive()
            toRev2 = PartialPlanarMechanics.FixedTranslation(r=[0.2,0.0])
            toMass1 = PartialPlanarMechanics.FixedTranslation(r=[0.1,0.0])
            mass1 = PartialPlanarMechanics.Body(m=m1,J=J1)
            rev2 = PartialPlanarMechanics.RevoluteJointActive()
            toTCP = PartialPlanarMechanics.FixedTranslation(r=[0.1,0.0])
            toMass2 = PartialPlanarMechanics.FixedTranslation(r=[0.05,0.0])
            mass2 = PartialPlanarMechanics.Body(m=m1,J=J1)
            tcp = TCP(eps=eps, mu_A= mu_A, mu_S=mu_S, vSlide=vSlide)
            flange1 = Rotational.Flange()
            flange2 = Rotational.Flange()
            zForceIn = Blocks.RealInput() 
            posXOut = Blocks.RealOutput()
            posYOut = Blocks.RealOutput()       
        end
        @equations begin
            connect(fixed.frame_a,rev1.frame_a)
            connect(rev1.frame_b,toRev2.frame_a)
            connect(rev1.frame_b,toMass1.frame_a)
            connect(toMass1.frame_b,mass1.frame_a)
            connect(toRev2.frame_b,rev2.frame_a)
            connect(rev2.frame_b,toTCP.frame_a)
            connect(rev2.frame_b,toMass2.frame_a)
            connect(toMass2.frame_b,mass2.frame_a)
            connect(toTCP.frame_b,tcp.frame_a)
            connect(zForceIn, tcp.zForceIn)
            connect(flange1,rev1.flange_a)
            connect(flange2,rev2.flange_a)
            connect(posXOut,tcp.posXOut)
            connect(posYOut,tcp.posYOut)
        end
    end

    @mtkmodel FrictionlessScaraMechanics begin
        @parameters begin
            m1 = 0.008, [description="mass of Link1"]
            J1 = 2.7e-5, [description="Inertia of Link1"]
            L1 = 0.2, [description="Length of Link1"]
            m2 = 0.004, [description="mass of Link2"]
            J2 = 3.5e-6, [description="Inertia of Link2"]
            L2 = 0.1, [description="Length of Link2"]
        end
        @components begin
            fixed = PartialPlanarMechanics.Fixed()
            rev1 = PartialPlanarMechanics.RevoluteJointActive()
            toRev2 = PartialPlanarMechanics.FixedTranslation(r=[0.2,0.0])
            toMass1 = PartialPlanarMechanics.FixedTranslation(r=[0.1,0.0])
            mass1 = PartialPlanarMechanics.Body(m=m1,J=J1)
            rev2 = PartialPlanarMechanics.RevoluteJointActive()
            toTCP = PartialPlanarMechanics.FixedTranslation(r=[0.1,0.0])
            toMass2 = PartialPlanarMechanics.FixedTranslation(r=[0.05,0.0])
            mass2 = PartialPlanarMechanics.Body(m=m1,J=J1)
            tcp = FrictionlessTCP()
            flange1 = Rotational.Flange()
            flange2 = Rotational.Flange()
            zForceIn = Blocks.RealInput() 
            posXOut = Blocks.RealOutput()
            posYOut = Blocks.RealOutput()       
        end
        @equations begin
            connect(fixed.frame_a,rev1.frame_a)
            connect(rev1.frame_b,toRev2.frame_a)
            connect(rev1.frame_b,toMass1.frame_a)
            connect(toMass1.frame_b,mass1.frame_a)
            connect(toRev2.frame_b,rev2.frame_a)
            connect(rev2.frame_b,toTCP.frame_a)
            connect(rev2.frame_b,toMass2.frame_a)
            connect(toMass2.frame_b,mass2.frame_a)
            connect(toTCP.frame_b,tcp.frame_a)
            connect(zForceIn, tcp.zForceIn)
            connect(flange1,rev1.flange_a)
            connect(flange2,rev2.flange_a)
            connect(posXOut,tcp.posXOut)
            connect(posYOut,tcp.posYOut)
        end
    end


    @mtkmodel ScaraMTKModel begin
        @structural_parameters begin
            inp1 = traj_1
            inp2 = traj_2
            inp3 = traj_3
        end

        @parameters begin
            m1 = 0.008, [description="mass of Link1"]
            J1 = 3.4e-4, [description="Inertia of Link1"]
            L1 = 0.2, [description="Length of Link1"]
            m2 = 0.004, [description="mass of Link2"]
            J2 = 4.35e-5, [description="Inertia of Link2"]
            L2 = 0.1, [description="Lengthof Link2"]
            elbow = -1.0, [description="elbow of Scara -1 - up, 1 down"]
            eps = 1e-6, [description="Stiction threshold"]
            vSlide = 0.001, [description = "Sliding velocity"]
            mu_A = 0.3, [description = "Friction coefficient at adhesion"]
            mu_S = 0.1, [description = "Friction coefficient at sliding"]
            R = 0.101, [description = "Motor Resistance [Ω]"]
            L = 2.66e-5, [description = "Motor Inductivity [H]"]
            k = 0.0115, [description = "Motor Constant [N.m/A]"]
            Jmotor = 0.0000119, [description = "Motor Inertia [kg.m2]"]
            v = 0.0, [description = "Viscous friction factor [Nm.s/rad]"]
            gearratio1 = 10.0, [description = "Gearratio for motor1"]
            gearratio2 = 5.0, [description = "Gearratio for motor2"]
            voltageLimit = 12.0, [description = "Limit of motor voltage input"]
            controllerGain2 = 4.423, [description = "Gain of P-Controller Motor1"]
            controllerGain1 = 7.1543, [description = "Gain of P-Controller Motor2"]
        end
        @components begin
            # input section
            inp1_traj = Blocks.TimeVaryingFunction(f = inp1)
            inp2_traj = Blocks.TimeVaryingFunction(f = inp2)
            inp3_traj = Blocks.TimeVaryingFunction(f = inp3)
            inv_kin = ScaraKinematic.ScaraInverseKinematics(L1=L1,L2=L2,elbow=elbow)
            # drive 1
            motor1 = BaseComponents.SignalDCMotor(R=R,L=L,k=k,J=Jmotor,v=v)
            gear1 = Rotational.IdealGear(ratio=gearratio1)
            angleSensor1 = Rotational.AngleSensor()
            feedback1 = Blocks.Feedback()
            gain1 = Blocks.Gain(k=controllerGain1)
            # drive 2
            motor2 = BaseComponents.SignalDCMotor(R=R,L=L,k=k,J=Jmotor,v=v)
            gear2 = Rotational.IdealGear(ratio = gearratio2)
            angleSensor2 = Rotational.AngleSensor()
            feedback2 = Blocks.Feedback()
            gain2 = Blocks.Gain(k=controllerGain2)
            scara = ScaraMechanics(m1=m1, J1=J1, L1=L1, m2=m2, J2=J2, L2=L2, eps=eps,
                vSlide=vSlide,mu_A=mu_A,mu_S=mu_S)
        end
        @equations begin
            connect(inp1_traj.output,inv_kin.inp1)
            connect(inp2_traj.output,inv_kin.inp2)
            # drive 1
            connect(inv_kin.out1,feedback1.input1)
            connect(feedback1.output, gain1.input)
            connect(gain1.output,motor1.vSignal)
            connect(motor1.flange_a,gear1.flange_a)
            connect(gear1.flange_b,angleSensor1.flange)
            connect(angleSensor1.phi,feedback1.input2)
            #drive 2
            connect(inv_kin.out2,feedback2.input1)
            connect(feedback2.output,gain2.input)
            connect(gain2.output, motor2.vSignal)
            connect(motor2.flange_a,gear2.flange_a)
            connect(gear2.flange_b,angleSensor2.flange)
            connect(angleSensor2.phi,feedback2.input2)
            # scara mechanics
            connect(gear1.flange_b,scara.flange1)
            connect(gear2.flange_b,scara.flange2)
            connect(inp3_traj.output,scara.zForceIn)
        end
    end


    @mtkmodel FrictionlessScaraMTKModel begin
        @structural_parameters begin
            inp1 = traj_1
            inp2 = traj_2
            inp3 = traj_3
        end

        @parameters begin
            m1 = 0.008, [description="mass of Link1"]
            J1 = 3.4e-4, [description="Inertia of Link1"]
            L1 = 0.2, [description="Length of Link1"]
            m2 = 0.004, [description="mass of Link2"]
            J2 = 4.35e-5, [description="Inertia of Link2"]
            L2 = 0.1, [description="Lengthof Link2"]
            elbow = -1.0, [description="elbow of Scara -1 - up, 1 down"]
            R = 0.101, [description = "Motor Resistance [Ω]"]
            L = 2.66e-5, [description = "Motor Inductivity [H]"]
            k = 0.0115, [description = "Motor Constant [N.m/A]"]
            Jmotor = 0.0000119, [description = "Motor Inertia [kg.m2]"]
            v = 0.0, [description = "Viscous friction factor [Nm.s/rad]"]
            gearratio1 = 10.0, [description = "Gearratio for motor1"]
            gearratio2 = 5.0, [description = "Gearratio for motor2"]
            voltageLimit = 12.0, [description = "Limit of motor voltage input"]
            controllerGain2 = 4.423, [description = "Gain of P-Controller Motor1"]
            controllerGain1 = 7.1543, [description = "Gain of P-Controller Motor2"]
        end
        @components begin
            # input section
            inp1_traj = Blocks.TimeVaryingFunction(f = inp1)
            inp2_traj = Blocks.TimeVaryingFunction(f = inp2)
            inp3_traj = Blocks.TimeVaryingFunction(f = inp3)
            inv_kin = ScaraKinematic.ScaraInverseKinematics(L1=L1,L2=L2,elbow=elbow)
            # drive 1
            motor1 = BaseComponents.SignalDCMotor(R=R,L=L,k=k,J=Jmotor,v=v)
            gear1 = Rotational.IdealGear(ratio=gearratio1)
            angleSensor1 = Rotational.AngleSensor()
            feedback1 = Blocks.Feedback()
            gain1 = Blocks.Gain(k=controllerGain1)
            # drive 2
            motor2 = BaseComponents.SignalDCMotor(R=R,L=L,k=k,J=Jmotor,v=v)
            gear2 = Rotational.IdealGear(ratio = gearratio2)
            angleSensor2 = Rotational.AngleSensor()
            feedback2 = Blocks.Feedback()
            gain2 = Blocks.Gain(k=controllerGain2)
            scara = FrictionlessScaraMechanics(m1=m1, J1=J1, L1=L1, m2=m2, J2=J2, L2=L2)
        end
        @equations begin
            connect(inp1_traj.output,inv_kin.inp1)
            connect(inp2_traj.output,inv_kin.inp2)
            # drive 1
            connect(inv_kin.out1,feedback1.input1)
            connect(feedback1.output, gain1.input)
            connect(gain1.output,motor1.vSignal)
            connect(motor1.flange_a,gear1.flange_a)
            connect(gear1.flange_b,angleSensor1.flange)
            connect(angleSensor1.phi,feedback1.input2)
            #drive 2
            connect(inv_kin.out2,feedback2.input1)
            connect(feedback2.output,gain2.input)
            connect(gain2.output, motor2.vSignal)
            connect(motor2.flange_a,gear2.flange_a)
            connect(gear2.flange_b,angleSensor2.flange)
            connect(angleSensor2.phi,feedback2.input2)
            # scara mechanics
            connect(gear1.flange_b,scara.flange1)
            connect(gear2.flange_b,scara.flange2)
            connect(inp3_traj.output,scara.zForceIn)
        end
    end



end