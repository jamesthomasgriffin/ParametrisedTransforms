// Dear Imgui includes taken from imgui_widgets.cpp
#ifndef IMGUI_DISABLE

// Used to add [] operator to ImVec4 (already present for ImVec2)
#define IM_VEC4_CLASS_EXTRA float& operator[](size_t ix) { return *(&(this->x) + ix); } \
float const& operator[](size_t ix) const { return *(&(this->x) + ix); }

#include "imgui.h"

#include <ctype.h>      // toupper
#if defined(_MSC_VER) && _MSC_VER <= 1500 // MSVC 2008 or earlier
#include <stddef.h>     // intptr_t
#else
#include <stdint.h>     // intptr_t
#endif

// Visual Studio warnings
#ifdef _MSC_VER
#pragma warning (disable: 4127)     // condition expression is constant
#pragma warning (disable: 4996)     // 'This function or variable may be unsafe': strcpy, strdup, sprintf, vsnprintf, sscanf, fopen
#if defined(_MSC_VER) && _MSC_VER >= 1922 // MSVC 2019 16.2 or later
#pragma warning (disable: 5054)     // operator '|': deprecated between enumerations of different types
#endif
#endif

// Clang/GCC warnings with -Weverything
#if defined(__clang__)
#if __has_warning("-Wunknown-warning-option")
#pragma clang diagnostic ignored "-Wunknown-warning-option"         // warning: unknown warning group 'xxx'                      // not all warnings are known by all Clang versions and they tend to be rename-happy.. so ignoring warnings triggers new warnings on some configuration. Great!
#endif
#pragma clang diagnostic ignored "-Wunknown-pragmas"                // warning: unknown warning group 'xxx'
#pragma clang diagnostic ignored "-Wold-style-cast"                 // warning: use of old-style cast                            // yes, they are more terse.
#pragma clang diagnostic ignored "-Wfloat-equal"                    // warning: comparing floating point with == or != is unsafe // storing and comparing against same constants (typically 0.0f) is ok.
#pragma clang diagnostic ignored "-Wformat-nonliteral"              // warning: format string is not a string literal            // passing non-literal to vsnformat(). yes, user passing incorrect format strings can crash the code.
#pragma clang diagnostic ignored "-Wsign-conversion"                // warning: implicit conversion changes signedness
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"  // warning: zero as null pointer constant                    // some standard header variations use #define NULL 0
#pragma clang diagnostic ignored "-Wdouble-promotion"               // warning: implicit conversion from 'float' to 'double' when passing argument to function  // using printf() is a misery with this as C++ va_arg ellipsis changes float to double.
#pragma clang diagnostic ignored "-Wenum-enum-conversion"           // warning: bitwise operation between different enumeration types ('XXXFlags_' and 'XXXFlagsPrivate_')
#pragma clang diagnostic ignored "-Wdeprecated-enum-enum-conversion"// warning: bitwise operation between different enumeration types ('XXXFlags_' and 'XXXFlagsPrivate_') is deprecated
#pragma clang diagnostic ignored "-Wimplicit-int-float-conversion"  // warning: implicit conversion from 'xxx' to 'float' may lose precision
#elif defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wpragmas"                  // warning: unknown option after '#pragma GCC diagnostic' kind
#pragma GCC diagnostic ignored "-Wformat-nonliteral"        // warning: format not a string literal, format string not checked
#pragma GCC diagnostic ignored "-Wclass-memaccess"          // [__GNUC__ >= 8] warning: 'memset/memcpy' clearing/writing an object of type 'xxxx' with no trivial copy-assignment; use assignment or value-initialization instead
#endif


#ifndef IMGUI_DEFINE_MATH_OPERATORS
#define IMGUI_DEFINE_MATH_OPERATORS
#endif

#include "imgui_internal.h"

#include "imgui_pt_widgets.h"

// Fill in the operations ImGui internal does not provide
static inline ImVec4 operator*(const ImVec4& lhs, const float rhs) { return ImVec4(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs, lhs.w * rhs); }
static inline ImVec4 operator/(const ImVec4& lhs, const float rhs) { return ImVec4(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs, lhs.w / rhs); }
static inline ImVec4 operator/(const ImVec4& lhs, const ImVec4& rhs) { return ImVec4(lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z, lhs.w / rhs.w); }
static inline ImVec4& operator*=(ImVec4& lhs, const float rhs) { lhs.x *= rhs; lhs.y *= rhs; lhs.z *= rhs; lhs.w *= rhs; return lhs; }
static inline ImVec4& operator/=(ImVec4& lhs, const float rhs) { lhs.x /= rhs; lhs.y /= rhs; lhs.z /= rhs; lhs.w /= rhs; return lhs; }
static inline ImVec4& operator+=(ImVec4& lhs, const ImVec4& rhs) { lhs.x += rhs.x; lhs.y += rhs.y; lhs.z += rhs.z; lhs.w += rhs.w; return lhs; }
static inline ImVec4& operator-=(ImVec4& lhs, const ImVec4& rhs) { lhs.x -= rhs.x; lhs.y -= rhs.y; lhs.z -= rhs.z; lhs.w -= rhs.w; return lhs; }
static inline ImVec4& operator*=(ImVec4& lhs, const ImVec4& rhs) { lhs.x *= rhs.x; lhs.y *= rhs.y; lhs.z *= rhs.z; lhs.w *= rhs.w; return lhs; }
static inline ImVec4& operator/=(ImVec4& lhs, const ImVec4& rhs) { lhs.x /= rhs.x; lhs.y /= rhs.y; lhs.z /= rhs.z; lhs.w /= rhs.w; return lhs; }
static inline float dot(ImVec2 const& a, ImVec2 const& b) { return a.x * b.x + a.y * b.y; }
static inline ImVec2 operator*(ImMat2 const& A, ImVec2 const& v) {
    return A.col1 * v.x + A.col2 * v.y;
}
static inline float dot(ImVec4 const& a, ImVec4 const& b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }
static inline ImVec4 operator*(ImMat4 const& A, ImVec4 const& v) {
    return A.col1 * v.x + A.col2 * v.y + A.col3 * v.z + A.col4 * v.w;
}


namespace ImGui
{
    /*
    This is a widget that generalises slider widgets from 1D to 2D/3D.

    In a classical slider widget there is a slider button whose position is
    linked to the value of a parameter.	Dragging the slider changes the value of
    the parameter.

    With the approach below the user defines control points which can depend on a
    number of different parameters, dragging these control points then changes
    the parameters.  The dependency between parameters and the position of a
    control point is fully customiseable via a stack of simple transforms that
    each depend on a parameter.	For example a transform may rotate by the angle
    specified by a given parameter.

    This is most easily understood by interacting with the examples.

    The position of the control point is specified by the stack of transformations
    which in turn depends on the parameters.  The tricky part is the reverse
    process, how to adjust the parameters to move the control point with the mouse
    input.  In general it wont be possible to match the control point position and
    mouse position exactly, but we can minimise the distance between them.

    This is achieved by calculating first and second derivatives of the control
    point position with respect to the parameter values, then applying the
    multi-variate Newton-Raphson method	to update the parameter values.

    The user of these classes does not have to worry about the implementation,
    they only need to specify the desired transformations and untransformed
    control point positions.
    */


    bool DraggableView::begin(const char* name, ImVec2 const& size, ImVec4 const& border_col)
    {
        IM_ASSERT(active == false);  // We shouldn't start our view twice

        ImGuiWindow* window = GetCurrentWindow();
        if (window->SkipItems)
            return false;

        m_size = size;
        m_cp_count = 0;

        ImGuiContext& g = *GImGui;
        //const ImGuiStyle& style = g.Style;
        const ImGuiID id = window->GetID(name);

        // Bounding box used for clipping 3D view and optional border
        m_initial_cursor_pos = GetCursorScreenPos();
        ImRect bb(window->DC.CursorPos, m_initial_cursor_pos + size);
        if (border_col.w > 0.0f)
            bb.Max += ImVec2(2, 2);
        ItemSize(bb);
        if (!ItemAdd(bb, 0))
            return false;

        m_final_cursor_pos = GetCursorScreenPos();

        if (border_col.w > 0.0f)
        {
            window->DrawList->AddRect(bb.Min, bb.Max, GetColorU32(border_col), 0.0f);
        }
        PushClipRect(bb.Min, bb.Max, true);

        active = true;

        return true;
    }

    void DraggableView::end()
    {
        if (!active)  // If we didn't start the view then there is nothing to do
            return;

        active = false;

        {   // Set the Dear ImGui state to where it should be after the view widget
            ImGui::SetCursorScreenPos(m_final_cursor_pos);

            PopClipRect();

            if (m_saved_control_point_id)
                SetActiveID(m_saved_control_point_id, GetCurrentWindow());
        }

        // Reset the saved control point, NB z-position does not need reseting
        m_saved_control_point_id = 0;

        // If there is a change in the deferral stack then apply it
        if (m_deferred_step_position >= 0)
            updateParameters(m_deferred_step);

        // Reset deferral stack state
        m_deferral_stack_size = 1;
        m_deferred_step_position = -1;
        m_deferred_step = {};
    }

    void DraggableView::pushDeferralSlot()
    {
        if (!active)
            return;
        ++m_deferral_stack_size;
    }

    void DraggableView::popDeferralSlot()
    {
        if (!active)
            return;

        IM_ASSERT(m_deferral_stack_size > 0);  // We shouldn't pop an empty stack

        // If there is a change at the top of the stack then apply it
        if (m_deferred_step_position + 1 == m_deferral_stack_size)
        {
            updateParameters(m_deferred_step);
            m_deferred_step = {};
            m_deferred_step_position = -1;
        }

        // Reduce the size
        --m_deferral_stack_size;
    }

    void DraggableView::drawDerivative(ImVec2 const& pos, ImVec2 const& d, ImVec4 const& col, float thickness) const
    {
        ImDrawList* draw_list = GetWindowDrawList();
        ImVec2 start_vector = viewCoordsToScreenCoords(pos);
        ImVec2 end_vector = viewCoordsToScreenCoords(pos + d);
        draw_list->AddLine(start_vector, end_vector, GetColorU32(col), thickness);
    }

    ImVec2 DraggableView::getMouseInViewCoords() const
    {
        ImVec2 mouse_pos = GetMousePos();
        return screenCoordsToViewCoords(mouse_pos);
    }

    ImVec2 DraggableView::viewCoordsToScreenCoords(ImVec2 const& v) const
    {
        ImVec2 result{ v };
        result += ImVec2(1, 1);
        result *= m_size / 2;
        result += m_initial_cursor_pos;
        return result;
    }

    ImVec2 DraggableView::screenCoordsToViewCoords(ImVec2 const& p) const
    {
        ImVec2 result{ p };
        result -= m_initial_cursor_pos;
        result /= m_size / 2;
        result -= ImVec2(1, 1);
        return result;
    }

    bool DraggableView::createControlPoint(ImVec2 const& marker_pos, ImGuiControlPointFlags const& flags, ImGuiButtonFlags const& button_flags, float marker_width, ImVec4 marker_col, float z_order) 
    {
        IM_ASSERT(active);  // We shouldn't have got this far with an inactive context

        // An Invisible button cannot have a size of 0, so limit it to 1.0f
        marker_width = (marker_width < 1.0f) ? 1.0f : marker_width;
        ImVec2 marker_size(marker_width, marker_width);

        ImRect marker_bb(marker_pos - marker_size / 2, marker_pos + marker_size / 2);

        m_last_control_point_position = marker_pos;
        m_last_control_point_radius = marker_width / 2;
        m_last_step = {};

        if (z_order < -1.0f)
        {
            // The control point is behind us so create a dummy icon we cannot interact with
            Dummy({ 0, 0 });
            return false;
        }

        bool cp_activated = false;            

        {   // Create invisible button at marker and return cursor to its original position
            char buf[32];
            sprintf(buf, "Marker %d", m_cp_count++);

            ImVec2 saved_cursor_screen_pos = GetCursorScreenPos();
            SetCursorScreenPos(marker_bb.Min);

            // Use an invisible button to handle hovering / clicking / dragging
            InvisibleButton(buf, marker_size, button_flags);

            SetItemAllowOverlap();

            if (IsItemActivated())  // Occurs only when first activated
            {
                cp_activated = true;
                if (flags & ImGuiControlPointFlags_Circular) {
                    // Perform test that mouse is over circle contained in bbox
                    ImVec2 vect_to_center = ImGui::GetMousePos() - (marker_bb.Min + marker_bb.Max) * 0.5f;
                    vect_to_center /= marker_width / 2;  // Divide by radius
                    float r2 = vect_to_center.x * vect_to_center.x + vect_to_center.y * vect_to_center.y;
                    if (r2 > 1.0f)
                    {
                        cp_activated = false;
                        ClearActiveID();
                            
                    }
                }
            }
            SetCursorScreenPos(saved_cursor_screen_pos);

            if (cp_activated)
            {
                saveControlPoint(z_order);
            }
        }

        // Draw the marker
        if ((flags & ImGuiControlPointFlags_DrawControlPointMarkers) && marker_col.w > 0.0f)
        {
            if (flags & ImGuiControlPointFlags_Circular) 
            {
                if (IsItemActive())
                    GetWindowDrawList()->AddCircleFilled((marker_bb.Min + marker_bb.Max) / 2, marker_width / 2, GetColorU32(marker_col));
                else
                    GetWindowDrawList()->AddCircle((marker_bb.Min + marker_bb.Max) / 2, marker_width / 2, GetColorU32(marker_col));
            } 
            else 
            {
                if (IsItemActive())
                    GetWindowDrawList()->AddRectFilled(marker_bb.Min, marker_bb.Max, GetColorU32(marker_col));
                else
                    GetWindowDrawList()->AddRect(marker_bb.Min, marker_bb.Max, GetColorU32(marker_col), 0.0f);
            }
        }
        return IsItemActive();
    }

    bool DraggableView::saveControlPoint(float z_order)
    {
        // If there is already a control point infront of this one then return
        if (m_saved_control_point_id && (z_order > m_saved_control_point_z_value))
            return false;

        m_saved_control_point_id = GetItemID();
        m_saved_control_point_z_value = z_order;
        return true;
    }

    void DraggableView::updateParameters(SavedParameterChanges const& changes)
    {
        updateParameters(changes, ImGuiControlPointFlags_ApplyParamChangesImmediately);
    }

    void DraggableView::updateParameters(SavedParameterChanges const& changes, ImGuiControlPointFlags const& flags)
    {
        m_last_step = changes;
        if (!(flags & ImGuiControlPointFlags_ApplyParamChangesImmediately)) {
            IM_ASSERT(m_deferral_stack_size > 0);  // There should be a slot in which to place our deferred changes
            m_deferred_step = changes;
            m_deferred_step_position = m_deferral_stack_size - 1;  // Step is at the end of the stack
        }
        else if (!(flags & ImGuiControlPointFlags_DoNotChangeParams)) {
            for (auto const& change : changes) {
                if (change.parameter)
                    *(change.parameter) += change.change;
            }
        }
    }



    ImVec2 Draggable2DView::apply(ImVec2 const& pos)
    {
        ImVec2 transformed_pos = m_transforms.applyToVector(pos);
        return viewCoordsToScreenCoords({ transformed_pos[0], transformed_pos[1] });
    }

    template<unsigned int D>
    void Draggable2DView::bringTogether(ImVec2 pos, ImVec2 p, std::array<float*, D> const& free_parameters, ImGuiControlPointFlags const& flags)
    {
        /*
        Minimises Q = |p - transform(pos)|^2 over the set of D parameters.
        When p is closest to Mv, the derivative with respect to a will be zero,
        so we use the Newton-Raphson method on the derivative to update the parameters.
        */

        // D is the number of parameters, there is nothing to do if D is zero
        if (D == 0)
            return;

        // Rather than directly apply the change from the Newton-Raphson method, first multiply
        // by this factor, slower convergence but smoother and less likely to jitter around
        constexpr float softening_factor = 0.5f;

        // A jet with given parameters, position and 0 first and second derivatives
        JetDeg2<ImVec2, float*, D> jet{ { free_parameters, pos } };

        DistanceSquaredToVector<float, ImVec2, float> Q{ p };
        auto output_jet = m_transforms.applyTo2Jet(jet);
        JetDeg2<float, float*, D> Qderivatives = Q.applyStaticallyTo2Jet(0, output_jet);

        std::array<float, D> steps = Qderivatives.applyConjugateGradientMethod();

        SavedParameterChanges changes{};
        for (int i = 0; i < D; ++i) {
            changes[i] = { free_parameters[i], -softening_factor * steps[i] };
        }
        updateParameters(changes, flags);
    }

    bool Draggable2DView::controlPoint(ImVec2 const& pos, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float importance)
    {
        return controlPoint<0>(pos, {}, flags, button_flags, marker_radius, marker_col, importance);
    }

    bool Draggable2DView::controlPoint(ImVec2 const& pos, float* free_param, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float importance)
    {
        return controlPoint<1>(pos, { free_param }, flags, button_flags, marker_radius, marker_col, importance);
    }

    bool Draggable2DView::controlPoint(ImVec2 const& pos, float* free_param1, float* free_param2, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float importance)
    {
        return controlPoint<2>(pos, { free_param1, free_param2 }, flags, button_flags, marker_radius, marker_col, importance);
    }

    template<unsigned int D>
    bool Draggable2DView::controlPoint(ImVec2 const& pos, std::array<float*, D> free_parameters, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float importance)
    {
        if (!active)  // We cannot add a control point to a view that hasn't been started
            return false;

        ImVec2 transformed_pos = m_transforms.applyToVector(pos);
        ImVec2 screen_coords = viewCoordsToScreenCoords({ transformed_pos[0], transformed_pos[1] });

        // Calculate marker size in pixels, using w coord to scale the size
        float marker_width = 2 * marker_radius;
        if (!(flags & ImGuiControlPointFlags_SizeInPixels))
            marker_width *= m_size.x;

        if (importance <= 0)
            importance = 0.0001f;

        // We want the control point to move the marker if it is the only one active (multiple can only be active at first activation)
        if (createControlPoint(screen_coords, flags, button_flags, marker_width, marker_col, 1.0 / importance) && IsItemActive() && !IsItemActivated())
        {
            ImVec2 mouse_pos = getMouseInViewCoords();

            bringTogether(pos, { mouse_pos[0], mouse_pos[1] }, free_parameters, flags);

            if (flags & ImGuiControlPointFlags_DrawParamDerivatives) {
                // Draw tangent vector(s)
                JetDeg1<ImVec2, float*, D> jet = m_transforms.applyTo1Jet<D>(JetDeg1<ImVec2, float*, D>{ free_parameters, pos });

                ImVec2 pos{ jet.position[0], jet.position[1] };
                for (auto const& d : jet.derivatives)
                    drawDerivative(pos, { d[0], d[1] }, marker_col);
            }

            return true;
        }
        return false;
    }

    void Draggable2DView::pushMatrix(ImMat2 const& M) { m_transforms.pushMatrix(M); }

    void Draggable2DView::pushRotation(float* parameter) { m_transforms.pushRotation(parameter); }
    void Draggable2DView::pushRotation(float angle) { m_transforms.pushRotation(angle); }

    void Draggable2DView::pushScale(float* parameter) { m_transforms.pushScale(parameter); }
    void Draggable2DView::pushScale(float factor) { m_transforms.pushScale(factor); }
    void Draggable2DView::pushScaleX(float* parameter) { m_transforms.pushScaleX(parameter); }
    void Draggable2DView::pushScaleX(float factor) { m_transforms.pushScaleX(factor); }
    void Draggable2DView::pushScaleY(float* parameter) { m_transforms.pushScaleY(parameter); }
    void Draggable2DView::pushScaleY(float factor) { m_transforms.pushScaleY(factor); }

    void Draggable2DView::pushTranslation(float* parameter, ImVec2 const& dir) { m_transforms.pushTranslation(parameter, dir); }
    void Draggable2DView::pushTranslation(float factor, ImVec2 const& dir) { m_transforms.pushTranslation(factor, dir); }
    void Draggable2DView::pushTranslationAlongX(float* parameter) { m_transforms.pushTranslationAlongX(parameter); }
    void Draggable2DView::pushTranslationAlongX(float factor) { m_transforms.pushTranslationAlongX(factor); }
    void Draggable2DView::pushTranslationAlongY(float* parameter) { m_transforms.pushTranslationAlongY(parameter); }
    void Draggable2DView::pushTranslationAlongY(float factor) { m_transforms.pushTranslationAlongY(factor); }


    ImVec2 Draggable3DView::apply(ImVec4 const& pos)
    {
        ImVec4 transformed_pos = m_transforms.applyToVector(pos);
        ImVec2 view_coords{ transformed_pos[0] / transformed_pos[3], transformed_pos[1] / transformed_pos[3] };
        return viewCoordsToScreenCoords(view_coords);
    }


    template<unsigned int D>
    void Draggable3DView::bringTogether(ImVec4 pos, ImVec4 mouse_in_view_coords, std::array<float*, D> const& free_parameters, ImGuiControlPointFlags const& flags)
    {
        // Minimises Q(a) = |p - Mv|^2 over a where Mv is the transformed vector, and
        // M(a) is the parametrised transformation with a the free parameter.
        // When p is closest to Mv, the derivative with respect to a will be
        // zero, so we use the Newton-Raphson method on the derivative.


        // Rather than directly apply the change from the Newton-Raphson method, first multiply
        // by this factor, slower convergence but smoother and less likely to jitter around
        constexpr float softening_factor = 0.5f;

        // A jet with given parameters, position and 0 first and second derivatives
        JetDeg2<ImVec4, float*, D> jet{ { free_parameters, pos } };

        // We need to project from 4 coords to 3 coords via p --> p / p.w
        ProjectionTransform<float, ImVec4, 3> projection{};
        DistanceSquaredToVector<float, ImVec4, float> Q{ mouse_in_view_coords };
        auto output_jet = projection.applyStaticallyTo2Jet(0, m_transforms.applyTo2Jet(jet));
        auto Qderivatives = Q.applyStaticallyTo2Jet(0, output_jet);

        std::array<float, D> steps = Qderivatives.applyConjugateGradientMethod();

        SavedParameterChanges changes{};
        for (int i = 0; i < D; ++i) {
            changes[i] = { free_parameters[i], -softening_factor * steps[i] };
        }
        updateParameters(changes, flags);
    }

    bool Draggable3DView::controlPoint(ImVec4 const& pos, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col)
    {
        return controlPoint<0>(pos, {}, flags, button_flags, marker_radius, marker_col);
    }

    bool Draggable3DView::controlPoint(ImVec4 const& pos, float* free_param, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col)
    {
        return controlPoint<1>(pos, { free_param }, flags, button_flags, marker_radius, marker_col);
    }

    bool Draggable3DView::controlPoint(ImVec4 const& pos, float* free_param1, float* free_param2, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col)
    {
        return controlPoint<2>(pos, { free_param1, free_param2 }, flags, button_flags, marker_radius, marker_col);
    }

    void Draggable3DView::pushMatrix(ImMat4 const& M) { m_transforms.pushMatrix(M); }

    void Draggable3DView::pushPerspectiveMatrix(float fov, float aspect_ratio, float z_near, float z_far)
    {
        float tan_half_fov = tanf(fov / 2);
        pushMatrix(ImMat4{
            {1 / (aspect_ratio * tan_half_fov), 0, 0, 0},
            {0, 1 / tan_half_fov, 0, 0},
            {0, 0, -(z_near + z_far) / (z_far - z_near), -1},
            {0, 0, -2 * z_far * z_near / (z_far - z_near), 0}
            });
    }

    void Draggable3DView::pushLookAtMatrix(Vec3<float> const& eye, Vec3<float> const& center, Vec3<float> const& up)
    {
        Vec3<float> const f{ normalize(center - eye) };
        Vec3<float> const s{ normalize(cross(f, up)) };
        Vec3<float> const u{ cross(s, f) };
        pushMatrix(ImMat4{
            { s.x, u.x, -f.x, 0 },
            { s.y, u.y, -f.y, 0 },
            { s.z, u.z, -f.z, 0 },
            {-dot(s, eye), -dot(u, eye), dot(f, eye), 1.0f } });
    }

    void Draggable3DView::pushRotationAboutAxis(float* parameter, Vec3<float> const& axis) { m_transforms.pushRotationAboutAxis(parameter, axis); }
    void Draggable3DView::pushRotationAboutAxis(float angle, Vec3<float> const& axis) { m_transforms.pushRotationAboutAxis(angle, axis); }
    void Draggable3DView::pushRotationAboutX(float* parameter) { m_transforms.pushRotationAboutX(parameter); }
    void Draggable3DView::pushRotationAboutX(float angle) { m_transforms.pushRotationAboutX(angle); }
    void Draggable3DView::pushRotationAboutY(float* parameter) { m_transforms.pushRotationAboutY(parameter); }
    void Draggable3DView::pushRotationAboutY(float angle) { m_transforms.pushRotationAboutY(angle); }
    void Draggable3DView::pushRotationAboutZ(float* parameter) { m_transforms.pushRotationAboutZ(parameter); }
    void Draggable3DView::pushRotationAboutZ(float angle) { m_transforms.pushRotationAboutZ(angle); }

    void Draggable3DView::pushScale(float* parameter) { m_transforms.pushScale(parameter); }
    void Draggable3DView::pushScale(float factor) { m_transforms.pushScale(factor); }
    void Draggable3DView::pushScaleX(float* parameter) { m_transforms.pushScaleX(parameter); }
    void Draggable3DView::pushScaleX(float factor) { m_transforms.pushScaleX(factor); }
    void Draggable3DView::pushScaleY(float* parameter) { m_transforms.pushScaleY(parameter); }
    void Draggable3DView::pushScaleY(float factor) { m_transforms.pushScaleY(factor); }
    void Draggable3DView::pushScaleZ(float* parameter) { m_transforms.pushScaleZ(parameter); }
    void Draggable3DView::pushScaleZ(float factor) { m_transforms.pushScaleZ(factor); }

    void Draggable3DView::pushTranslation(float* parameter, ImVec4 const& dir) { m_transforms.pushTranslation(parameter, dir); }
    void Draggable3DView::pushTranslation(float factor, ImVec4 const& dir) { m_transforms.pushTranslation(factor, dir); }
    void Draggable3DView::pushTranslationAlongX(float* parameter) { m_transforms.pushTranslationAlongX(parameter); }
    void Draggable3DView::pushTranslationAlongX(float factor) { m_transforms.pushTranslationAlongX(factor); }
    void Draggable3DView::pushTranslationAlongY(float* parameter) { m_transforms.pushTranslationAlongY(parameter); }
    void Draggable3DView::pushTranslationAlongY(float factor) { m_transforms.pushTranslationAlongY(factor); }
    void Draggable3DView::pushTranslationAlongZ(float* parameter) { m_transforms.pushTranslationAlongZ(parameter); }
    void Draggable3DView::pushTranslationAlongZ(float factor) { m_transforms.pushTranslationAlongZ(factor); }

    template<unsigned int D>
    bool Draggable3DView::controlPoint(ImVec4 const& pos, std::array<float*, D> free_parameters, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col)
    {
        if (!active)  // We cannot add a control point to a view that hasn't been started
            return false;

        ImVec4 transformed_pos = m_transforms.applyToVector(pos);
        ImVec2 screen_coords = viewCoordsToScreenCoords(
            ImVec2{ transformed_pos[0] / transformed_pos[3], transformed_pos[1] / transformed_pos[3] }
        );

        float z = transformed_pos[2] / transformed_pos[3];

        // Calculate marker size in pixels, using w coord to scale the size
        float marker_width = 2 * marker_radius;
        if(!(flags & ImGuiControlPointFlags_SizeInPixels))
            marker_width *= m_size.x;
        if(!(flags & ImGuiControlPointFlags_FixedSize))
            marker_width /= (transformed_pos.w);

        if (createControlPoint(screen_coords, flags, button_flags, marker_width, marker_col, z)) {
            if ((D > 0) && IsItemActive() && !IsItemActivated()) 
            {
                ImVec2 mouse_pos = getMouseInViewCoords();
                ImVec4 mouse_in_view_coords{ mouse_pos[0], mouse_pos[1], transformed_pos.z / transformed_pos.w, 1.0f };
                
                bringTogether(pos, mouse_in_view_coords, free_parameters, flags);

                if (flags & ImGuiControlPointFlags_DrawParamDerivatives) 
                {
                    // Draw tangent vector(s)
                    ProjectionTransform<float, ImVec4, 3> projection{};
                    JetDeg1<ImVec4, float*, D> jet = projection.applyStaticallyTo1Jet<D>(0,
                        m_transforms.applyTo1Jet<D>(JetDeg1<ImVec4, float*, D>{ free_parameters, pos })
                        );
                    ImVec2 pos{ jet.position[0], jet.position[1] };
                    for (auto const& d : jet.derivatives)
                        drawDerivative(pos, { d[0], d[1] }, marker_col);
                }
                return true;
            }
        }
        return false;
    }

    Draggable3DView g_current_context{};

    bool BeginControlPointView(const char* name, ImVec2 const& size, ImVec4 const& border_color)
    {
        return g_current_context.begin(name, size, border_color);
    }

    void EndControlPointView()
    {
        g_current_context.end();
    }

    void PushDeferralSlot()
    {
        g_current_context.pushDeferralSlot();
    }

    void PopDeferralSlot()
    {
        g_current_context.popDeferralSlot();
    }

    ImVec2 GetLastControlPointPosition() {
        return g_current_context.getLastControlPointPosition();
    }

    float GetLastControlPointRadius() {
        return g_current_context.getLastControlPointRadius();
    }

    bool ControlPoint(ImVec2 const& pos, float* free_parameter, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float z_order)
    {
        return g_current_context.controlPoint({ pos.x, pos.y, z_order, 1 }, free_parameter, flags, button_flags, marker_radius, marker_col);
    }
    bool ControlPoint(ImVec2 const& pos, float* free_parameter1, float* free_parameter2, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float z_order)
    {
        return g_current_context.controlPoint({ pos.x, pos.y, z_order, 1 }, free_parameter1, free_parameter2, flags, button_flags, marker_radius, marker_col);
    }
    bool ControlPoint(ImVec2 const& pos, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float z_order)
    {
        return g_current_context.controlPoint({pos.x, pos.y, z_order, 1}, flags, button_flags, marker_radius, marker_col);
    }

    bool ControlPoint(ImVec4 const& pos, float* free_parameter, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col)
    {
        return g_current_context.controlPoint(pos, free_parameter, flags, button_flags, marker_radius, marker_col);
    }
    bool ControlPoint(ImVec4 const& pos, float* free_parameter1, float* free_parameter2, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col)
    {
        return g_current_context.controlPoint(pos, free_parameter1, free_parameter2, flags, button_flags, marker_radius, marker_col);
    }
    bool ControlPoint(ImVec4 const& pos, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col)
    {
        return g_current_context.controlPoint(pos, flags, button_flags, marker_radius, marker_col);
    }

    void PushMatrix(ImMat4 const& M) { g_current_context.pushMatrix(M); }
    void PushPerspectiveMatrix(float fov, float aspect_ratio, float z_near, float z_far) { g_current_context.pushPerspectiveMatrix(fov, aspect_ratio, z_near, z_far); }
    void PushLookAtMatrix(Vec3<float> const& eye, Vec3<float> const& center, Vec3<float> const& up) { g_current_context.pushLookAtMatrix(eye, center, up); }
    void PushMatrix(ImMat2 const& M) { 
        ImMat4 N {
            { M.col1.x, M.col1.y, 0, 0 },
            { M.col2.x, M.col2.y, 0, 0 },
            { 0, 0, 1, 0 },
            { 0, 0, 0, 1 } };
        g_current_context.pushMatrix(N); 
    }

    void PushRotationAboutAxis(float* parameter, Vec3<float> const& axis) { g_current_context.pushRotationAboutAxis(parameter, axis); }
    void PushRotationAboutAxis(float angle, Vec3<float> const& axis) { g_current_context.pushRotationAboutAxis(angle, axis); }
    void PushRotationAboutX(float* parameter) { g_current_context.pushRotationAboutX(parameter); }
    void PushRotationAboutX(float angle) { g_current_context.pushRotationAboutX(angle); }
    void PushRotationAboutY(float* parameter) { g_current_context.pushRotationAboutY(parameter); }
    void PushRotationAboutY(float angle) { g_current_context.pushRotationAboutY(angle); }
    void PushRotationAboutZ(float* parameter) { g_current_context.pushRotationAboutZ(parameter); }
    void PushRotationAboutZ(float angle) { g_current_context.pushRotationAboutZ(angle); }
    void PushRotation(float* parameter) { PushRotationAboutZ(parameter); }
    void PushRotation(float angle) { PushRotationAboutZ(angle); }

    void PushScale(float* parameter) { g_current_context.pushScale(parameter); }
    void PushScale(float factor) { g_current_context.pushScale(factor); }
    void PushScaleX(float* parameter) { g_current_context.pushScaleX(parameter); }
    void PushScaleX(float factor) { g_current_context.pushScaleX(factor); }
    void PushScaleY(float* parameter) { g_current_context.pushScaleY(parameter); }
    void PushScaleY(float factor) { g_current_context.pushScaleY(factor); }
    void PushScaleZ(float* parameter) { g_current_context.pushScaleZ(parameter); }
    void PushScaleZ(float factor) { g_current_context.pushScaleZ(factor); }

    void PushTranslation(float* parameter, ImVec4 const& v) { g_current_context.pushTranslation(parameter, v); }
    void PushTranslation(float value, ImVec4 const& v) { g_current_context.pushTranslation(value, v); }
    void PushTranslation(float* parameter, ImVec2 const& v) { g_current_context.pushTranslation(parameter, { v.x, v.y, 0, 0 }); }
    void PushTranslation(float value, ImVec2 const& v) { g_current_context.pushTranslation(value, { v.x, v.y, 0, 0 }); }
    void PushTranslationAlongX(float* parameter) { g_current_context.pushTranslationAlongX(parameter); }
    void PushTranslationAlongX(float value) { g_current_context.pushTranslationAlongX(value); }
    void PushTranslationAlongY(float* parameter) { g_current_context.pushTranslationAlongY(parameter); }
    void PushTranslationAlongY(float value) { g_current_context.pushTranslationAlongY(value); }
    void PushTranslationAlongZ(float* parameter) { g_current_context.pushTranslationAlongZ(parameter); }
    void PushTranslationAlongZ(float value) { g_current_context.pushTranslationAlongZ(value); }

    void PopTransform() { g_current_context.popTransform(); }
    void ClearTransforms() { g_current_context.clearTransforms(); }


    
    void TranslationWidget(ImMat4* M, float axes_size, ImGuiControlPointFlags control_point_flags, bool include_planes)
    {
        ImVec4 t{ 0, 0, 0, 0 };  // The translation vector (only three coords are used)
        ImGui::PushTranslationAlongX(&t.x);
        ImGui::PushTranslationAlongY(&t.y);
        ImGui::PushTranslationAlongZ(&t.z);

        // We don't change any parameters in the control points, instead wait until the deferral slot is popped
        ImGui::PushDeferralSlot();

        ImGui::ControlPoint(M->col4, control_point_flags, 0, 0.025f, { 1, 1, 1, 1 });
        ImGui::ControlPoint(*M * ImVec4{ axes_size, 0, 0, 1 }, &t.x, control_point_flags, 0, 0.025f, { 1, 0, 0, 1 });
        ImGui::ControlPoint(*M * ImVec4{ 0, axes_size, 0, 1 }, &t.y, control_point_flags, 0, 0.025f, { 0, 1, 0, 1 });
        ImGui::ControlPoint(*M * ImVec4{ 0, 0, axes_size, 1 }, &t.z, control_point_flags, 0, 0.025f, { 0, 0, 1, 1 });

        if (include_planes) {
            ImGui::ControlPoint(*M * ImVec4{ axes_size / 2, axes_size / 2, 0, 1 }, &t.x, &t.y, control_point_flags, 0, 0.02f, { 1, 1, 0, 1 });
            ImGui::ControlPoint(*M * ImVec4{ 0, axes_size / 2, axes_size / 2, 1 }, &t.y, &t.z, control_point_flags, 0, 0.02f, { 0, 1, 1, 1 });
            ImGui::ControlPoint(*M * ImVec4{ axes_size / 2, 0, axes_size / 2, 1 }, &t.z, &t.x, control_point_flags, 0, 0.02f, { 1, 0, 1, 1 });
        }

        ImGui::PopTransform();  // Z-translation
        ImGui::PopTransform();  // Y-translation
        ImGui::PopTransform();  // X-translation

        // This will apply any changes made to t as intended, without the unintended consequence of
        // applying any previously made changes.
        ImGui::PopDeferralSlot();

        // Update the matrix if any of our parameters are non-zero
        M->col4 += *M * t;
    }

}  // End ImGui namespace

#endif  // IMGUI_DISABLE