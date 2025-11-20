#
# Copyright (c) 2025 Andreas Hofmann
# Licensed under the MIT license. See LICENSE file in the project root for details.
#

# Basic components required for ScaraRobot example
# partially based on the Modelica library Servomechanisms https://github.com/afrhu/Servomechanisms

module ScaraKinematic

    using ModelingToolkit
    using ModelingToolkit: t_nounits as t, D_nounits as D
    using ModelingToolkitStandardLibrary.Blocks
    using NaNMath

    @mtkmodel ScaraInverseKinematics begin
        @parameters begin
            L1 = 1.0, [description = "Length of first Link"]
            L2 = 1.0, [description = "Length of second Link"]
            elbow = 1.0, [description = "elbow direction: -1 elbow up, 1 elbow down"]
        end
        @variables begin
            u1(t)
            u2(t)
            y1(t)
            y2(t)
        end
        @components begin
            inp1 = Blocks.RealInput()
            inp2 = Blocks.RealInput()
            out1 = Blocks.RealOutput()
            out2 = Blocks.RealOutput()
        end
        @equations begin
            u1 ~ inp1.u
            u2 ~ inp2.u
            y1 ~ out1.u
            y2 ~ out2.u
            y2 ~ atan(elbow * NaNMath.sqrt(1 - ((u1^2 + u2^2 - L1^2 - L2^2) / (2 * L1 * L2))^2), (u1^2 + u2^2 - L1^2 - L2^2) / (2 * L1 * L2))
            y1 ~ atan(u2, u1) - atan(L2 * sin(y2), L1 + L2 * cos(y2));
        end
    end
    
end