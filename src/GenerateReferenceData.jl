#
# Copyright (c) 2025 Andreas Hofmann
# Licensed under the MIT license. See LICENSE file in the project root for details.
#

module GenerateReferenceData

    export generateReferenceData

    using ..ScaraRobot
    using ModelingToolkit
    using DataInterpolations
    using DelimitedFiles
    using DifferentialEquations
    using HDF5

    
    function generateReferenceData(inputdata::String, simname::String, h5file::String;samplingrate=1e-4)

        function getScaraInputs(path_to_data)
            trajs = readdlm(path_to_data)
            
            interpolation_objects = []
            for i in 2:4
                push!(interpolation_objects,LinearInterpolation(trajs[:,i],trajs[:,1]; extrapolation = ExtrapolationType.Linear))
            end
            return interpolation_objects

        end

        # training data
        (inp1,inp2,inp3) = getScaraInputs(inputdata )

        @mtkcompile scara = ScaraRobot.ScaraMTKModel(inp1=inp1, inp2=inp2, inp3=inp3)

        # result quantities
        state_quantities = unknowns(scara)
        
        input_quantities = [
            (inp1,"inp1"),
            (inp2,"inp2"),
            (inp3,"inp3")
        ]

        output_quantities = [
            (scara.scara.tcp.vx,"TCP.vx"),
            (scara.scara.tcp.vy,"TCP.vy"),
            (scara.scara.tcp.x,"TCP.x"),
            (scara.scara.tcp.y,"TCP.y"),
            (scara.scara.tcp.fx,"TCP.FrictionForceX"),
            (scara.scara.tcp.fy,"TCP.FrictionForceY")
        ]

        

        tspan = (inp1.t[1],inp1.t[end])

        u0 = [
            scara.scara.rev1.phi => 0.0,
            scara.scara.rev2.phi => 0.0,
            scara.scara.rev1.w => 0.0,
            scara.scara.rev2.w => 0.0,
            scara.motor1.inductor.i => 0.0,
            scara.motor2.inductor.i => 0.0
        ]

        prob = ODEProblem(scara,u0,tspan)

        sol = solve(prob,Rodas5P(),saveat=samplingrate, tstops=inp1.t)

        h5open(h5file,"w") do file
            gr = create_group(file,simname)

            for (i,entry) in enumerate(state_quantities)
                item_group = create_group(gr,string(entry))
                item_group["value"] = sol(sol.t,Val{0};idxs=i).u
                item_group["time"] = sol.t
                der_group = create_group(gr,"der_"*string(entry))
                value = sol(sol.t,Val{1};idxs=i).u
                der_group["value"] = value
                der_group["time"] = sol.t
            end

            for entry in input_quantities
                item_group = create_group(gr,entry[2])
                item_group["value"] = entry[1].u
                item_group["time"] = entry[1].t
            end

            for entry in output_quantities
                item_group = create_group(gr,entry[2])
                item_group["value"] = sol[entry[1]]
                item_group["time"] = sol.t
                
            end

        end

        return string.(state_quantities)

    end

end