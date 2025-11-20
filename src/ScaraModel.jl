#
# Copyright (c) 2025 Andreas Hofmann
# Licensed under the MIT license. See LICENSE file in the project root for details.
#

module ScaraModel

    #using ModelingToolkit

    export ScaraRobot
    export TrainableScaraRobot

    include("BaseComponents.jl")
    include("PartialPlanarMechanics.jl")
    include("ScaraKinematic.jl")
    include("SymbolicInputFunctions.jl")
    include("ScaraRobot.jl")
    include("TrainableScaraRobot.jl")
    include("GenerateReferenceData.jl")

end # module ScaraModel