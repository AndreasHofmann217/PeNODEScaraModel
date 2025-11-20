#
# Copyright (c) 2025 Andreas Hofmann
# Licensed under the MIT license. See LICENSE file in the project root for details.
#


module SymbolicInputFunctions

    using ModelingToolkit

    export traj_1, traj_2, traj_3


    # the ScaraRobot has 3 inputs: joint1, joint2 and force into the plane
    # since the inputs will be replaced later, we define some simple symbolic functions here, that will later be manually replaced after code generation
    function traj_1(t)
        return 1.0
    end
    function traj_2(t)
        return 2.0
    end
    function traj_3(t)
        return 3.0
    end
    @register_symbolic traj_1(t)
    @register_symbolic traj_2(t)
    @register_symbolic traj_3(t)


end