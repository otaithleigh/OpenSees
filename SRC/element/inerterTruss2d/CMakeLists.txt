#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        clutchInerterTruss2d.cpp
        inerterTruss2d.cpp
    PUBLIC
        clutchInerterTruss2d.h
        inerterTruss2d.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
