# [Reference frames creation and transformations](@id tutorial_06_frames)


```julia
using Basic, ReferenceFrameRotations
```

## Introduction

`Basic` provides the capability to create, expand and efficienty differentiate Frame Systems. 
A frame system allows you to compute the relative position and orientation of a registered object
with respect to any other registered objects at any time. This includes time derivatives of motion
up to order 3 (jerk) and well integrated with `ForwardDiff`, allowing to compute partial 
derivatives of motion with respect to _any_ parameter.

`Frames` is a rather generic framework which can be extended and configured by the user. The Frame 
System, indeed, is the _ultimate goal_ of the `Basic` computational environment. Therefore, 
the aim of this tutorial is to provide an overview of the core functionalities of a Frame System 
together with some relevant informations relative to its structure and conceptual modeling, 
in such a way any user would be capable to understand, use and extend it in the proper way.

## The `FrameSystem`

The entry-point in the `Frames` computational environment, is the `FrameSystem` object. A 
`FrameSystem` manages a collection of points and axes in the form of `FramePointNode`s and 
`FrameAxesNode`s, respectively. These two kind of _nodes_, as you can see from the 
[Points Graphs Tutorial](@ref tutorial_04_points) and the [Axes Graphs Tutorial](@ref tutorial_05_axes),
are organized in a graph form and represent two precise entities:

- **Axes**: defines an orientation in space. These are related each other by means of a `Rotation` 
    transformation which relate one axes to a parent axes in a certain time interval.
- **Points**: defines a location in space. These are related each other by means of a `Translation`
    transformation which relate one point to a parent point in a particular axes in a certain 
    time interval.

Any node can have several transformations defining their orientation or position with 
respect to other `FrameAxesNode` or `FramePointNode`, each applicable during a particular 
time period. Moreover, as you have seen in the dedicated tutorials:

- Nodes can be created independently of each other (by means of `@axes` or `@point` macros).
- They shall be **registered** within the `FrameSystem` to be used, connecting with a parent by
    means of a transformation. This step can be usually performed by means of a dedicated method,
    i.e. `add_axes_xxxx!` or `add_point_yyyy!`, where `xxxx` and `yyyy` shall be specified 
    depending on the particular axes or point type.

!!! note
    Don't confuse the concept of frame of reference and axes! In fact, while in general these
    two concept are mixed together, they have a precise meaning within `Basic`:
    
    - Axes define coordinate systems, which are a mathematical concept;
    - A frame of reference is a physical concept related to the state of motion of something.

    In practice, we use a _coordinate system_ to specify a frame of reference because it is 
    convenient from a computational perspective but to each frame of reference there should 
    in theory correspond a family of an infinite numbers of coordinate systems comoving 
    together. 
    
    To clarify more, the velocity of a point with respect to a physical three-dimensional 
    body is an absolute concept (in Newtonian mechanics) but to express it as a set of 3 
    components we need to choose a coordinate system. A natural choice is to use a 
    coordinate system comoving with that body but it is also possible to use a coordinate 
    system rotating with respect to that body.

    It is normal to express quantities computed in a particular frame of reference in a 
    different coordinate system. For example, one might solve the motion of spacecraft in 
    an inertial frame of reference centered on the Earth, but express its velocity coordinates 
    along axes rotating with the Earth. This is not the same as computing the velocity 
    of the spacecraft in the Earth rotating system.





