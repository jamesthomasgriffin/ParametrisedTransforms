# NB

I thought this was neat enough to develop further, so have started from scratch with a different model for the transformation stack and a different focus, only on matrix transformations.  The new code is far more efficient and will be available soon, it will be called ImControl.

# ParametrisedTransforms

A small template header library for parametrised transforms, with examples for Dear ImGui.

![arm demo](https://user-images.githubusercontent.com/2971239/103889425-51301700-50de-11eb-8e90-e20fde6e0885.gif)

Read the source files for documentation.  To run the demos start with your chosen Dear ImGui implementation example and call each of the three declared functions immediately before or after the ShowDemo function call.

If you have any questions about the code then please contact me or raise an issue.  Furthermore if you have any suggestions for transforms you want included then again raise an issue.  Contributions to the code are of course welcome.

You may notice some strange behaviour when dragging certain control points around.  This is probably a result of the underlying geometry, the simplest way to mitigate such issues is to restrict the domain of the parameters.

### Tests

There are some tests I used when developing the transforms, they check that the transforms do what I expected them to do and they perform a numerical check of the hard-coded formal derivatives.  To run them you need to link the glm library and the googletest library.
