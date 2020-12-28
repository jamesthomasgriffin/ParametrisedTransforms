#pragma once

#include <array>
#include <memory>
#include <utility>
#include <vector>

// Next steps
// * move glm's rotation function into this code
// * add more examples
// * split the transformations and the contexts into separate files
// * add code that applies Newton-Raphson to a jet

/*

Parametrised Transforms

Our transforms are a function that takes a vector and returns another vector, for
example a rotation by a quarter turn.  We intend that the transforms are twice
differentiable, meaning they have first and second derivatives.  They do not have
to be affine like a matrix transform, twice differentiability is the only restriction.
Another example is projection of a 4-vector: sending p to p / p.w.

A simple parametrised transform has in addition a single parameter, for example
a rotation by t radians has a parameter t.

A general parametrised transform may have a list of parameters, an example with 
parameters (r, s, t) might be to rotate around the x-axis by r radians, then to 
translate in the y-axis by a distance s, before rotating in the z-axis by t radians.
Our parametrised transforms will all be expressed as a composition of a number
of simple transforms.


Discussion of Usage

The code relies on runtime polymorphism and so to transform a vector we have to
reference a long list of virtual methods, applying each in turn.  Unlike matrix
transforms the list of transformations can not be collapsed to a single matrix and
then applied.  The intended usage is to apply a variety of transforms, decided on 
at runtime, to a small number of points.

If performance is an issue and one does need to apply transforms to many points 
there is a simple trick that allows this.  Just write a vector class that holds
many vectors in parallel and such that v[0] is the array of x-values and so on.
In particular it is not necessary for v[i] to be of the type V::value_type, it can
itself be a vector.  This can allow a transformation to act on many points at once,
with all the performance optimisations, simd etc. that entails.


Linear Algebra Template Parameters

Each of the templated classes requires one or two vector types, denoted by V and W.
They also require a scalar type, scalar_t.  The classes must also implement, 

* scalar multiplication, 
* vector addition, 
* initialisation by a single scalar, a copy constructor, copy assignment,
* access to coordinates via [] operator

The other type of class we need for some of our templates is denoted P, the only 
requirement of this class is that it has the comparison operator ==.  It is used 
as a label for parameters, often it will be convenient that P be a pointer to a
scalar and dereferencing the pointer will give the scalar value of the parameter.

We use a matrix class denoted M_t in a few of the templates, these must be equipped
with a left multiplication with the vector class V giving values in W where 
appropriate.

The remaining template parameters are unsigned integers, D denotes the number
of parameters, while I and J denote particular coordinates.
*/

/*

Jets

Jets are used in differential geometry, they represent a point along with
additional information in the tangent space at that point.  For example an
object of type JetDeg1<V, P, D> uses a single vector of type V to represent 
position, and D vectors of type V to represent D tangent vectors pointing out from
that point, each of these tangent vectors is labelled by an object of type P.

The idea is that the tangent vector corresponding to the label p is the derivative of
the position with respect to a change in the parameter labelled by p.

An object of type JetDeg2<V, P, D> is an object of type JetDeg1<V, P, D> along with
D choose 2 additional vectors, this time labelled by pairs of objects of type P, these
are the second derivatives.

These jets are convenient because they contain all information required to be
transformed by the transformation, so just as the transformation transforms points,
it also transforms jets, the transformation on jets is specified by the chain rule.
*/
template<typename V, typename P, unsigned int D>
struct JetDeg1 {
	// The variables differentiated against, not their values but a label to compare against,
	// pointers to variables can work well. 
	std::array<P, D> variables;

	// The position of the jet
	V position{};

	// The derivative with respect to each of the variables
	std::array<V, D> derivatives{};

	// Not needed but put in to align with second derivative accessors below
	V& derivative(int i) { return derivatives[i]; }
	V const& derivative(int i) const { return derivatives[i]; }

	JetDeg1 projective_jet(int W) const {
		// First transform the position
		auto pw = position[W];
		JetDeg1 result{ variables, position / pw };

		// Now for the derivatives, we apply the quotient rule
		for (int i = 0; i < D; ++i) {
			result.derivatives[i] = (derivatives[i] * pw - position * derivatives[i][W]) / (pw * pw);
		}
		return result;
	}
};

// A 2-jet is a 1-jet plus information about 2nd derivatives
template<typename V, typename P, unsigned int D>
struct JetDeg2 : public JetDeg1<V, P, D>
{
	// The second derivative with respect to each pair of variables
	std::array<V, (D* (D + 1)) / 2> second_derivatives{};

	V& second_derivative(int i, int j) {
		if (i > j)
			std::swap(i, j);
		return second_derivatives[i + (j * (j + 1)) / 2];
	}
	V const& second_derivative(int i, int j) const {
		if (i > j)
			std::swap(i, j);
		return second_derivatives[i + (j * (j + 1)) / 2];
	}

	std::array<V, D> applyConjugateGradientMethod() const {
		// Compute first_derivative / second_derivative using the conjugate gradient method.
		// Typically this will be applied when V is a scalar type.

		constexpr float error_tolerance{ 0.01f };

		// Since we start with x = 0, the initialisation may look slightly different to standard
		std::array<V, D> x{};  // Start with x = 0
		std::array<V, D> r{ derivatives };
		std::array<V, D> d{ derivatives };

		V delta_new{};
		for (auto const& v : r)
			delta_new += v * v;

		V const delta_0 = delta_new;

		for (int i = 0; (i < D) && (delta_new > error_tolerance * error_tolerance * delta_0); ++i)
		{
			std::array<V, D> q{};
			applySecondDerivative(q, d);

			V alpha{};
			for (int j = 0; j < D; ++j)
				alpha += d[j] * q[j];
			alpha = delta_new / alpha;

			for (int j = 0; j < D; ++j) {
				x[j] += alpha * d[j];
				r[j] -= alpha * q[j];
			}

			V delta_old{ delta_new };

			delta_new = V{};
			for (auto const& v : r)
				delta_new += v * v;

			V beta = delta_new / delta_old;

			for (int j = 0; j < D; ++j)
				d[j] = r[j] + beta * d[j];
		}
		return x;
	}

	void applySecondDerivative(std::array<V, D>& out, std::array<V, D> const& x, V const& f = 1) const {
		// out <== out + f second_derivative * x
		for (int i = 0; i < D; ++i) {
			for (int j = 0; j < D; ++j) {
				out[i] += second_derivative(i, j) * x[j];
			}
		}
	}
};

// An abstract base class for a 1 parameter family of transformations p --> T(t, p) acting 
// on vectors p of type V and producing a vector of type W.
// The vector valued function (t, p) --> T(t, p) must be twice differentiable.
template<typename S, typename V, typename W>
class BaseSimpleTransform
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = W;

	virtual ~BaseSimpleTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const = 0;

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const = 0;
	virtual out_vector_t D(scalar_t t, in_vector_t const& p, in_vector_t const& v) const = 0;

	virtual out_vector_t D2(scalar_t t, in_vector_t const& p, in_vector_t const& v, in_vector_t const& w) const = 0;
	virtual out_vector_t Ddt(scalar_t t, in_vector_t const& p, in_vector_t const& v) const = 0;
	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const = 0;

	template<unsigned int D, typename P>
	JetDeg1<out_vector_t, P, D> applyTo1Jet(scalar_t t, P variable, JetDeg1<in_vector_t, P, D> const& J) const {
		// Initialise with the variables from the input and the position
		JetDeg1<out_vector_t, P, D> result{ 
			J.variables, 
			this->applyToVector(t, J.position) 
		};

		// Now fill in the derivative
		for (int i = 0; i < D; ++i)
		{
			// Apply the static transformation to the derivative
			result.derivative(i) = this->D(t, J.position, J.derivative(i));

			// If the variable matches, there is an additional term from the parameter derivative
			if (variable == J.variables[i])
				result.derivative(i) += dt(t, J.position);
		}
		return result;
	}

	template<unsigned int D, typename P>
	JetDeg2<out_vector_t, P, D> applyTo2Jet(scalar_t t, P variable, JetDeg2<in_vector_t, P, D> const& J) const {
		// Initialise using the application to the base 1-jet
		JetDeg2<out_vector_t, P, D> result{
			this->applyTo1Jet<D, P>(t, variable, J)
		};

		// Now fill in the second derivatives
		for (int i = 0; i < D; ++i)	{
			for (int j = i; j < D; ++j) {
				result.second_derivative(i, j) = this->D(t, J.position, J.second_derivative(i, j));
				result.second_derivative(i, j) += this->D2(t, J.position, J.derivative(i), J.derivative(j));
				if (variable == J.variables[i])
				{
					result.second_derivative(i, j) += this->Ddt(t, J.position, J.derivative(j));
				}
				if (variable == J.variables[j])
				{
					result.second_derivative(i, j) += this->Ddt(t, J.position, J.derivative(i));
				}
				if ((variable == J.variables[i]) && (variable == J.variables[j]))
				{
					result.second_derivative(i, j) += this->dt2(t, J.position);
				}
			}
		}
		return result;
	}

	// Used if the parameter t and parameters P are known to be independent
	template<unsigned int D, typename P>
	JetDeg1<out_vector_t, P, D> applyStaticallyTo1Jet(scalar_t t, JetDeg1<in_vector_t, P, D> const& J) const {
		// Initialise with the variables from the input and the position
		JetDeg1<out_vector_t, P, D> result{
			J.variables,
			this->applyToVector(t, J.position)
		};

		// Now fill in the derivative
		for (int i = 0; i < D; ++i)
		{
			// Apply the static transformation to the derivative
			result.derivative(i) = this->D(t, J.position, J.derivative(i));
		}
		return result;
	}

	// Used if the parameter t and parameters P are known to be independent
	template<unsigned int D, typename P>
	JetDeg2<out_vector_t, P, D> applyStaticallyTo2Jet(scalar_t t, JetDeg2<in_vector_t, P, D> const& J) const {
		// Initialise using the application to the base 1-jet
		JetDeg2<out_vector_t, P, D> result{
			this->applyStaticallyTo1Jet<D, P>(t, J)
		};

		// Now fill in the second derivatives
		for (int i = 0; i < D; ++i) {
			for (int j = i; j < D; ++j) {
				result.second_derivative(i, j) = this->D(t, J.position, J.second_derivative(i, j));
				result.second_derivative(i, j) += this->D2(t, J.position, J.derivative(i), J.derivative(j));
			}
		}
		return result;
	}
};

// If T is a complete descendent of BaseSimpleTransform, then this class is an abstract base class for 
// transformations of the form T(g(t)).  It's descendents only need implement g and its first and second
// derivatives.
template<typename T>
class ReparametrisationOfSimpleTransform : public T
{
public:
	ReparametrisationOfSimpleTransform(T const& transform) : T{ transform } {}

	// scalar_t and vector_t are inherited from T
	using scalar_t = typename T::scalar_t;
	using in_vector_t = typename T::in_vector_t;
	using out_vector_t = typename T::out_vector_t;

	// Define the reparametrisation and its first and second derivatives
	virtual scalar_t g(scalar_t t) const = 0;
	virtual scalar_t dgdt(scalar_t t) const = 0;
	virtual scalar_t d2gdt2(scalar_t t) const = 0;

	// Functions required to implement a simple transform are constructed from those of T and the above functions
	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) { return T::applyToVector(g(t), p); }

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) { return T::dt(g(t), p) * dgdt(t); }
	virtual out_vector_t D(scalar_t t, in_vector_t const& p, in_vector_t const& v) { return T::D(g(t), p, v); }

	virtual out_vector_t D2(scalar_t t, in_vector_t const& p, in_vector_t const& v, in_vector_t const& w) { return T::D2(g(t), p, v, w);  }
	virtual out_vector_t Ddt(scalar_t t, in_vector_t const& p, in_vector_t const& v) { return T::Ddt(g(t), p, v) * dgdt(t); }
	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) { return T::dt2(g(t), p) * dgdt(t) * dgdt(t) + T::dt(g(t), p) * d2gdt2(t); }
};

template<typename T>
class LinearReparametrisation : public ReparametrisationOfSimpleTransform<T>
{
public:
	using scalar_t = typename T::scalar_t;

	LinearReparametrisation(scalar_t const& a, scalar_t const& b, T const& transform) :
		ReparametrisationOfSimpleTransform<T>{ transform }, m_a{ a }, m_b{ b } {}
	LinearReparametrisation(scalar_t const& a, T const& transform) :
		ReparametrisationOfSimpleTransform<T>{ transform }, m_a{ a }, m_b{ 0 } {}
	virtual ~LinearReparametrisation() {}

	virtual scalar_t g(scalar_t t) const { return m_a * t + m_b; }
	virtual scalar_t dgdt(scalar_t t) const { return m_a; }
	virtual scalar_t d2gdt2(scalar_t t) const { return 0; }

private:
	scalar_t m_a;
	scalar_t m_b;
};

// A class that represents a family of transformations with multiple parameters.
// These transformations are expressed as a product of simple transformations.
template<typename S, typename V, typename P>
struct ParametrisedTransform
{
	using scalar_t = S;
	using vector_t = V;
	using transformation_t = BaseSimpleTransform<scalar_t, vector_t, vector_t>;

	class SimpleTransform {
	public:
		template<typename T>
		SimpleTransform(T trans, P param) : m_transformation{ std::make_shared<T>(trans) }, m_parameter{ param }, m_value{}, m_should_dereference{ true } {}

		template<typename T>
		SimpleTransform(T trans, scalar_t val) : m_transformation{ std::make_shared<T>(trans) }, m_parameter{}, m_value{ val }, m_should_dereference{ false } {}

		scalar_t getValue() const {	return m_should_dereference ? *m_parameter : m_value; }
		P const& getParameter() const { return m_parameter; }
		transformation_t const& getTransformation() const { return *m_transformation; }

	private:
		std::shared_ptr<transformation_t> m_transformation;
		P m_parameter;
		scalar_t m_value;
		bool m_should_dereference{ true };
	};

	std::vector<SimpleTransform> simple_transforms{};

	vector_t applyToVector(vector_t const& p) const {
		vector_t result{ p };
		for (auto it = simple_transforms.crbegin();	it != simple_transforms.crend(); ++it) {  // Iterate backwards through transforms
			result = it->getTransformation().applyToVector(it->getValue(), result);
		}
		return result;
	}

	template<unsigned int D>
	JetDeg1<vector_t, P, D> applyTo1Jet(JetDeg1<vector_t, P, D> const& p) const {
		JetDeg1<vector_t, P, D> result{ p }; 
		for (auto it = simple_transforms.crbegin();	it != simple_transforms.crend(); ++it) {
			result = it->getTransformation().applyTo1Jet<D>(it->getValue(), it->getParameter(), result);
		}
		return result;
	}

	template<unsigned int D>
	JetDeg2<vector_t, P, D> applyTo2Jet(JetDeg2<vector_t, P, D> const& p) const {
		JetDeg2<vector_t, P, D> result{ p };
		for (auto it = simple_transforms.crbegin(); it != simple_transforms.crend(); ++it) {
			result = it->getTransformation().applyTo2Jet<D>(it->getValue(), it->getParameter(), result);
		}
		return result;
	}

	ParametrisedTransform operator*(ParametrisedTransform const& b) const {
		ParametrisedTransform output{ *this };  // Perform copy using default copy constructor, Qu would a deep copy be more appropriate?
		output.simple_transforms.insert(output.simple_transforms.end(), b.simple_transforms.begin(), b.simple_transforms.end());
		return output;
	}
	ParametrisedTransform& operator*=(ParametrisedTransform const& b) {
		simple_transforms.insert(simple_transforms.end(), b.simple_transforms.begin(), b.simple_transforms.end());
		return *this;
	}


	template<unsigned int D>
	JetDeg1<vector_t, P, D> dt(vector_t const& p, std::array<P, D> const& parameters) const {
		// Apply the transform to the jet at input vector p, with given parameters and zero derivatives
		return this->applyToVector(
			JetDeg1<vector_t, P, D>{ parameters, p }
		);
	}

	vector_t D(vector_t const& p, vector_t const& v) const {
		vector_t output{ v };
		vector_t q{ p };
		for (auto const& transform : simple_transforms) {
			output = transform.getTransformation().D(q, output);
			q = transform.getTransformation().applyToVector(q);
		}
		return output;
	}

	vector_t D2(vector_t const& p, vector_t const& v, vector_t const& w) const;
	template<unsigned int D>
	vector_t Ddt(vector_t const& p, vector_t const& v, std::array<P, D> const& parameters) const;
	template<unsigned int D>
	JetDeg2<vector_t, P, D> dt2(vector_t const& p, std::array<P, D> const& parameters) const {
		// Apply the transform to the jet at input vector p, with given parameters and zero first and second derivatives
		return this->applyToVector(
			JetDeg2<vector_t, P, D>{ parameters, p }
		);
	}
};

template<typename S, typename V, typename W>
class StaticTransform : public BaseSimpleTransform<S, V, W>
{
public:
	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const { return out_vector_t{}; }
	virtual out_vector_t Ddt(scalar_t t, in_vector_t const& p, in_vector_t const& v) const { return out_vector_t{}; }
	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const { return out_vector_t{}; }
};

template<typename S, typename V, unsigned int J>
class ProjectionTransform : public StaticTransform<S, V, V>
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = V;

	virtual ~ProjectionTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const {
		return p / p[J];
	}

	virtual out_vector_t D(scalar_t t, in_vector_t const& p, in_vector_t const& v) const {
		return (v * p[J] - p * v[J]) / (p[J] * p[J]);
	}

	virtual out_vector_t D2(scalar_t t, in_vector_t const& p, in_vector_t const& v, in_vector_t const& w) const {
		return (p * (2 * v[J] * w[J]) - w * (p[J] * v[J]) - v * (p[J] * w[J])) / (p[J] * p[J] * p[J]);
	}
};

// For affine transforms, i.e. transforms acting on vectors via v --> M(t)v + m(t),
// the second spatial derivatives are zero.
// This allows for the following definition.
template<typename S, typename V, typename W>
class AffineFamilyOfSimpleTransforms : public BaseSimpleTransform<S, V, W>
{
public:
	//using in_vector_t = V;
	//using out_vector_t = W;

	~AffineFamilyOfSimpleTransforms() {};

	virtual out_vector_t D2(scalar_t t, in_vector_t const& p, in_vector_t const& v, in_vector_t const& w) const { return out_vector_t{}; }
};

// For linear transforms, i.e. transforms acting on vectors via v --> M(t)v, the spatial 
// derivatives are just the matrix M(t), i.e. the action on the spatial derivatives is the
// same as the action of the transform.
// This allows for the following definitions.
template<typename S, typename V, typename W>
class MatrixFamilyOfSimpleTransforms : public AffineFamilyOfSimpleTransforms<S, V, W>
{
public:
	~MatrixFamilyOfSimpleTransforms() {};
	virtual out_vector_t D(scalar_t t, in_vector_t const& p, in_vector_t const& v) const { return applyToVector(t, v); }
	virtual out_vector_t Ddt(scalar_t t, in_vector_t const& p, in_vector_t const& v) const { return dt(t, v); }
};

// A transform given by v --> Mv + m, with no dependence on a scalar parameter
template<typename S, typename V, typename W, typename M_t>
class StaticAffineTransform : public AffineFamilyOfSimpleTransforms<S, V, W>
{
public:
	StaticAffineTransform(M_t const& M, W const& m) : m_M{ M }, m_m{ m } {}
	explicit StaticAffineTransform(M_t const& M) : m_M{ M }, m_m{} {}
	~StaticAffineTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const { return m_M * p + m_m; }

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const { return out_vector_t{}; }
	virtual out_vector_t D(scalar_t t, in_vector_t const& p, in_vector_t const& v) const { return m_M * v; }

	virtual out_vector_t Ddt(scalar_t t, in_vector_t const& p, in_vector_t const& v) const { return out_vector_t{}; }
	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const { return out_vector_t{}; }

private:
	M_t m_M;  // Linear term
	W m_m;    // Constant term
};

// A transform given by matrix multiplication v -- > Mv with no dependence on a scalar parameter
template<typename S, typename V, typename W, typename M_t>
class StaticMatrixTransform : public MatrixFamilyOfSimpleTransforms<S, V, W>
{
public:
	explicit StaticMatrixTransform(M_t const& M) : m_M{ M } {}
	~StaticMatrixTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const { return m_M * p; }

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const { return out_vector_t{}; }

	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const { return out_vector_t{}; }

private:
	M_t m_M;  // Linear term
};

// A transform that acts by v --> v + t dir for a given direction vector, dir
template<typename S, typename V>
class TranslateDirectionTransform : public AffineFamilyOfSimpleTransforms<S, V, V>
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = V;

	TranslateDirectionTransform(out_vector_t const& dir) : m_dir{ dir } {}
	~TranslateDirectionTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const { return p + m_dir * t; }

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const { return m_dir; }

	virtual out_vector_t D(scalar_t t, in_vector_t const& p, in_vector_t const& v) const { return out_vector_t{ v }; }
	virtual out_vector_t Ddt(scalar_t t, in_vector_t const& p, in_vector_t const& v) const { return out_vector_t{}; }

	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const { return out_vector_t{}; }

private:
	out_vector_t m_dir;
};

// A transform that acts by scalar multiplication by the parameter, v --> t v
template<typename S, typename V>
class ScaleTransform : public MatrixFamilyOfSimpleTransforms<S, V, V>
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = V;

	~ScaleTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const {
		return p * t;
	}

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const {
		return p;
	}

	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const { return out_vector_t{}; }
};

// Acts via scaling but for projective geometry, i.e. it leaves a given coordinate unchanged
template<typename S, typename V, unsigned int J>
class ProjectiveScaleTransform : public MatrixFamilyOfSimpleTransforms<S, V, V>
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = V;
	~ProjectiveScaleTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const {
		out_vector_t result{ p * t };
		result[J] = p[J];
		return result;
	}

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const {
		out_vector_t result{ p };
		result[J] = 0;
		return result;
	}

	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const { return out_vector_t{}; }
};

// Translates by an axis vector, v --> v + t e_i for one of the axis vectors e_i
template<typename S, typename V, unsigned int axis>
class TranslateAlongAxisTransform : public AffineFamilyOfSimpleTransforms<S, V, V>
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = V;

	~TranslateAlongAxisTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const { 
		out_vector_t result{ p };
		result[axis] += t;
		return result;
	}

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const { 
		out_vector_t result{};
		result[axis] = 1;
		return result;
	}

	virtual out_vector_t D(scalar_t t, in_vector_t const& p, in_vector_t const& v) const { return out_vector_t{ v }; }
	virtual out_vector_t Ddt(scalar_t t, in_vector_t const& p, in_vector_t const& v) const { return out_vector_t{}; }

	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const { return out_vector_t{}; }
};

// Scales but along only one coordiate 
template<typename S, typename V, unsigned int axis>
class ScaleAlongAxisTransform : public MatrixFamilyOfSimpleTransforms<S, V, V>
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = V;

	~ScaleAlongAxisTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const {
		out_vector_t result{ p };
		result[axis] *= t;
		return result;
	}

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const {
		out_vector_t result{};
		result[axis] = p[axis];
		return result; 
	}

	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const { return out_vector_t{}; }
};

// Is specialised to rotations about a particular axis in the 3D case
template<typename S, typename V, unsigned int I, unsigned int J>
class RotateIJTransform : public MatrixFamilyOfSimpleTransforms<S, V, V>
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = V;

	~RotateIJTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const {
		out_vector_t result{ p };
		float c = cos(t), s = sin(t);
		result[I] = p[I] * c - p[J] * s;
		result[J] = p[I] * s + p[J] * c;
		return result;
	}

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const {
		out_vector_t result{};
		float c = cos(t), s = sin(t);
		result[I] = -p[I] * s - p[J] * c;
		result[J] = p[I] * c - p[J] * s;
		return result;
	}

	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const {
		out_vector_t result{};
		float c = cos(t), s = sin(t);
		result[I] = -p[I] * c + p[J] * s;
		result[J] = -p[I] * s - p[J] * c;
		return result;
	}
};

template<typename scalar_t, typename V>
using RotateXTransform = RotateIJTransform<scalar_t, V, 1, 2>;
template<typename scalar_t, typename V>
using RotateYTransform = RotateIJTransform<scalar_t, V, 2, 0>;
template<typename scalar_t, typename V>
using RotateZTransform = RotateIJTransform<scalar_t, V, 0, 1>;

// A skew transform, i.e. action of one of the basis matrices, v --> v + t E_ij v
// This is a translation in projective geometry
template<typename S, typename V, unsigned int I, unsigned int J>
class SkewTransform : public MatrixFamilyOfSimpleTransforms<S, V, V>
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = V;

	~SkewTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const {
		out_vector_t result{ p };
		result[I] += t * p[J];
		return result;
	}

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const {
		out_vector_t result{};
		result[I] = p[J];
		return result;
	}

	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const {
		return out_vector_t{};
	}
};

// See the definition for what this does, it's main use is in projective
// geometry, where for a vec3 v, SkewDirectionTransform<vec4, 3>{ {v, 0} }
// is a translation in the direction v.
template<typename S, typename V, unsigned int J>
class SkewDirectionTransform : public MatrixFamilyOfSimpleTransforms<S, V, V>
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = V;

	explicit SkewDirectionTransform(V const& v) : m_dir{ v } {};
	~SkewDirectionTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const {
		out_vector_t result{ p };
		result += m_dir * (t * p[J]);
		return result;
	}

	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const {
		return m_dir * p[J];
	}

	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const {
		return out_vector_t{};
	}

private:
	V m_dir;
};


// Returns the distance squared to a given point
template<typename S, typename V, typename W>
class DistanceSquaredToVector : public StaticTransform<S, V, W>
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = W;

	explicit DistanceSquaredToVector(V const& v) : m_vector{ v } {};
	~DistanceSquaredToVector() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const {
		in_vector_t diff{ p - m_vector };
		return dot(diff, diff) / 2;
	}

	virtual out_vector_t D(scalar_t t, in_vector_t const& p, in_vector_t const& v) const {
		return dot(v, p - m_vector);
	}

	virtual out_vector_t D2(scalar_t t, in_vector_t const& p, in_vector_t const& v, in_vector_t const& w) const {
		return dot(v, w);
	}

private:
	V m_vector;
};

// A very simple vector class used for rotations (which rely on 3 dimensions)
template<typename S>
struct Vec3 {
	using scalar_t = S;
	scalar_t x{}, y{}, z{};
	scalar_t& operator[](size_t ix) { return *(&x + ix); }
	scalar_t const& operator[](size_t ix) const { return *(&x + ix); }

	Vec3 operator+(Vec3 const& other) const { return Vec3{ x + other.x, y + other.y, z + other.z }; }
	Vec3 operator-(Vec3 const& other) const { return Vec3{ x - other.x, y - other.y, z - other.z }; }
	Vec3 operator*(float a) const { return Vec3{ a * x, a * y, a * z }; }
	Vec3 operator/(float a) const { return Vec3{ x / a, y / a, z / a }; }

	Vec3& operator+=(Vec3 const& other) { x += other.x; y += other.y; z += other.z; return *this; }
	Vec3& operator-=(Vec3 const& other) { x -= other.x; y -= other.y; z -= other.z; return *this; }
	Vec3& operator*=(float a) { x *= a; y *= a; z *= a; return *this; }
	Vec3& operator/=(float a) { x /= a; y /= a; z /= a; return *this; }
};

template<typename S> static inline S dot(Vec3<S> const& a, Vec3<S> const& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
template<typename S> static inline S length(Vec3<S> const& a) { return sqrt(dot(a, a)); }
template<typename S> static inline Vec3<S> normalize(Vec3<S> const& v) { return v / length(v); }
template<typename S> static inline Vec3<S> cross(Vec3<S> const& a, Vec3<S> const& b) { return Vec3<S>{ a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }


// Rotates around the given vector by 'angle' radians following the convention set in the glm library.
template<typename S, typename V>
class RotateAboutAxisTransform : public MatrixFamilyOfSimpleTransforms<S, V, V>
{
public:
	using scalar_t = S;
	using in_vector_t = V;
	using out_vector_t = V;

	RotateAboutAxisTransform(Vec3<scalar_t> const& axis) : m_axis{ normalize(axis) } {};
	~RotateAboutAxisTransform() {};

	virtual out_vector_t applyToVector(scalar_t t, in_vector_t const& p) const {
		// Adapted from the glm implementation of rotate
		scalar_t const c = cos(t);
		scalar_t const s = sin(t);

		// Calculate result in 3D
		Vec3<scalar_t> temp = m_axis * (1 - c);

		scalar_t rotate[3][3];  // Write the below as an initialiser list?
		rotate[0][0] = c + temp[0] * m_axis[0];  // Following column-order only because glm does
		rotate[0][1] = temp[0] * m_axis[1] + s * m_axis[2];
		rotate[0][2] = temp[0] * m_axis[2] - s * m_axis[1];

		rotate[1][0] = temp[1] * m_axis[0] - s * m_axis[2];
		rotate[1][1] = c + temp[1] * m_axis[1];
		rotate[1][2] = temp[1] * m_axis[2] + s * m_axis[0];

		rotate[2][0] = temp[2] * m_axis[0] + s * m_axis[1];
		rotate[2][1] = temp[2] * m_axis[1] - s * m_axis[0];
		rotate[2][2] = c + temp[2] * m_axis[2];

		out_vector_t result{ p };
		for (int i = 0; i < 3; ++i) {
			result[i] = {};
			for (int j = 0; j < 3; ++j)
				result[i] += rotate[j][i] * p[j];
		}
		return result;
	}
	virtual out_vector_t dt(scalar_t t, in_vector_t const& p) const {
		// Differentiated applyToVector
		scalar_t const c = cos(t);
		scalar_t const s = sin(t);

		Vec3<scalar_t> dtemp = m_axis * s;

		scalar_t dRotate[3][3];
		dRotate[0][0] = -s + dtemp[0] * m_axis[0];
		dRotate[0][1] = dtemp[0] * m_axis[1] + c * m_axis[2];
		dRotate[0][2] = dtemp[0] * m_axis[2] - c * m_axis[1];

		dRotate[1][0] = dtemp[1] * m_axis[0] - c * m_axis[2];
		dRotate[1][1] = -s + dtemp[1] * m_axis[1];
		dRotate[1][2] = dtemp[1] * m_axis[2] + c * m_axis[0];

		dRotate[2][0] = dtemp[2] * m_axis[0] + c * m_axis[1];
		dRotate[2][1] = dtemp[2] * m_axis[1] - c * m_axis[0];
		dRotate[2][2] = -s + dtemp[2] * m_axis[2];

		out_vector_t result{};
		for (int i = 0; i < 3; ++i) {
			result[i] = {};
			for (int j = 0; j < 3; ++j)
				result[i] += dRotate[j][i] * p[j];
		}
		return result;
	}
	virtual out_vector_t dt2(scalar_t t, in_vector_t const& p) const {
		// Second derivative of applyToVector
		scalar_t const a = t;
		scalar_t const c = cos(a);
		scalar_t const s = sin(a);

		Vec3<scalar_t> d2temp = m_axis * c;

		scalar_t d2Rotate[3][3];
		d2Rotate[0][0] = -c + d2temp[0] * m_axis[0];
		d2Rotate[0][1] = d2temp[0] * m_axis[1] - s * m_axis[2];
		d2Rotate[0][2] = d2temp[0] * m_axis[2] + s * m_axis[1];

		d2Rotate[1][0] = d2temp[1] * m_axis[0] + s * m_axis[2];
		d2Rotate[1][1] = -c + d2temp[1] * m_axis[1];
		d2Rotate[1][2] = d2temp[1] * m_axis[2] - s * m_axis[0];

		d2Rotate[2][0] = d2temp[2] * m_axis[0] - s * m_axis[1];
		d2Rotate[2][1] = d2temp[2] * m_axis[1] + s * m_axis[0];
		d2Rotate[2][2] = -c + d2temp[2] * m_axis[2];

		out_vector_t result{};
		for (int i = 0; i < 3; ++i) {
			result[i] = {};
			for (int j = 0; j < 3; ++j)
				result[i] += d2Rotate[j][i] * p[j];
		}
		return result;
	}

private:
	Vec3<scalar_t> const m_axis;
};

// Simple enhancement of ParameterisedTransform to add convenience functions
template<typename S, typename V, typename M_t, typename P>
struct ParametrisedTransformContext : public ParametrisedTransform<S, V, P>
{
	using scalar_t = S;
	using vector_t = V;
	using matrix_t = M_t;

	using SimpleTransform = typename ParametrisedTransform<scalar_t, vector_t, P>::SimpleTransform;

	template<typename Trans>
	void pushSingleParameterTransformation(P parameter) {
		simple_transforms.push_back(
			SimpleTransform(Trans{}, parameter)
		);
	}

	template<typename Trans, typename arg_type>
	void pushSingleParameterTransformationWithArgument(P parameter, arg_type const& argument) {
		simple_transforms.push_back(
			SimpleTransform(Trans{ argument }, parameter)
		);
	}

	template<typename Trans>
	void pushStaticTransformation(scalar_t value) {
		simple_transforms.push_back(
			SimpleTransform(Trans{}, value)
		);
	}

	template<typename Trans, typename arg_type>
	void pushStaticTransformationWithArgument(scalar_t value, arg_type const& argument) {
		simple_transforms.push_back(
			SimpleTransform(Trans{ argument }, value)
		);
	}

	void pushMatrix(matrix_t const& M) { pushStaticTransformationWithArgument<StaticMatrixTransform<scalar_t, vector_t, vector_t, matrix_t>, matrix_t>(0, M); }

	void pushRotationAboutAxis(P parameter, Vec3<scalar_t> const& axis) { pushSingleParameterTransformationWithArgument<RotateAboutAxisTransform<scalar_t, vector_t>, Vec3<scalar_t>>(parameter, axis); }
	void pushRotationAboutAxis(scalar_t angle, Vec3<scalar_t> const& axis) { pushStaticTransformationWithArgument<RotateAboutAxisTransform<scalar_t, vector_t>, Vec3<scalar_t>>(angle, axis); }
	void pushRotationAboutX(P parameter) { pushSingleParameterTransformation<RotateIJTransform<scalar_t, vector_t, 1, 2>>(parameter); }
	void pushRotationAboutX(scalar_t angle) { pushStaticTransformation<RotateIJTransform<scalar_t, vector_t, 1, 2>>(angle); }
	void pushRotationAboutY(P parameter) { pushSingleParameterTransformation<RotateIJTransform<scalar_t, vector_t, 2, 0>>(parameter); }
	void pushRotationAboutY(scalar_t angle) { pushStaticTransformation<RotateIJTransform<scalar_t, vector_t, 2, 0>>(angle); }
	void pushRotationAboutZ(P parameter) { pushSingleParameterTransformation<RotateIJTransform<scalar_t, vector_t, 0, 1>>(parameter); }
	void pushRotationAboutZ(scalar_t angle) { pushStaticTransformation<RotateIJTransform<scalar_t, vector_t, 0, 1>>(angle); }

	void pushScale(P parameter) { pushSingleParameterTransformation<ScaleTransform<scalar_t, vector_t>>(parameter); }
	void pushScale(scalar_t factor) { pushStaticTransformation<ScaleTransform<scalar_t, vector_t>>(factor); }
	void pushScaleX(P parameter) { pushSingleParameterTransformation<ScaleAlongAxisTransform<scalar_t, vector_t, 0>>(parameter); }
	void pushScaleX(scalar_t factor) { pushStaticTransformation<ScaleAlongAxisTransform<scalar_t, vector_t, 0>>(factor); }
	void pushScaleY(P parameter) { pushSingleParameterTransformation<ScaleAlongAxisTransform<scalar_t, vector_t, 1>>(parameter); }
	void pushScaleY(scalar_t factor) { pushStaticTransformation<ScaleAlongAxisTransform<scalar_t, vector_t, 1>>(factor); }
	void pushScaleZ(P parameter) { pushSingleParameterTransformation<ScaleAlongAxisTransform<scalar_t, vector_t, 2>>(parameter); }
	void pushScaleZ(scalar_t factor) { pushStaticTransformation<ScaleAlongAxisTransform<scalar_t, vector_t, 2>>(factor); }

	void pushTranslation(P parameter, vector_t v) { pushSingleParameterTransformationWithArgument<TranslateDirectionTransform<scalar_t, vector_t>, vector_t>(parameter, v); }
	void pushTranslation(scalar_t value, vector_t v) { pushStaticTransformationWithArgument<TranslateDirectionTransform<scalar_t, vector_t>, vector_t>(value, v); }
	void pushTranslationAlongX(P parameter) { pushSingleParameterTransformation<TranslateAlongAxisTransform<scalar_t, vector_t, 0>>(parameter); }
	void pushTranslationAlongX(scalar_t value) { pushStaticTransformation<TranslateAlongAxisTransform<scalar_t, vector_t, 0>>(value); }
	void pushTranslationAlongY(P parameter) { pushSingleParameterTransformation<TranslateAlongAxisTransform<scalar_t, vector_t, 1>>(parameter); }
	void pushTranslationAlongY(scalar_t value) { pushStaticTransformation<TranslateAlongAxisTransform<scalar_t, vector_t, 1>>(value); }
	void pushTranslationAlongZ(P parameter) { pushSingleParameterTransformation<TranslateAlongAxisTransform<scalar_t, vector_t, 2>>(parameter); }
	void pushTranslationAlongZ(scalar_t value) { pushStaticTransformation<TranslateAlongAxisTransform<scalar_t, vector_t, 2>>(value); }

	void popTransform() { simple_transforms.pop_back(); }

	void clearTransforms() { simple_transforms.clear(); }
};

template<typename S, typename V, typename M>
struct ParametrisedProjective3dTransformContext : public ParametrisedTransformContext<S, V, M, S*>
{
	using scalar_t = S;
	using vector_t = V;
	using P = scalar_t*;
	using matrix_t = M;

	void pushScale(P parameter) { pushSingleParameterTransformation<ScaleTransform<scalar_t, vector_t>>(parameter); }
	void pushScale(scalar_t factor) { pushStaticTransformation<ScaleTransform<scalar_t, vector_t>>(factor); }

	void pushTranslation(P parameter, vector_t v) { pushSingleParameterTransformationWithArgument<SkewDirectionTransform<scalar_t, vector_t, 3>, vector_t>(parameter, v); }
	void pushTranslation(scalar_t value, vector_t v) { pushStaticTransformationWithArgument<SkewDirectionTransform<scalar_t, vector_t, 3>, vector_t>(value, v); }
	void pushTranslationAlongX(P parameter) { pushSingleParameterTransformation<SkewTransform<scalar_t, vector_t, 0, 3>>(parameter); }
	void pushTranslationAlongX(scalar_t value) { pushStaticTransformation<SkewTransform<scalar_t, vector_t, 0, 3>>(value); }
	void pushTranslationAlongY(P parameter) { pushSingleParameterTransformation<SkewTransform<scalar_t, vector_t, 1, 3>>(parameter); }
	void pushTranslationAlongY(scalar_t value) { pushStaticTransformation<SkewTransform<scalar_t, vector_t, 1, 3>>(value); }
	void pushTranslationAlongZ(P parameter) { pushSingleParameterTransformation<SkewTransform<scalar_t, vector_t, 2, 3>>(parameter); }
	void pushTranslationAlongZ(scalar_t value) { pushStaticTransformation<SkewTransform<scalar_t, vector_t, 2, 3>>(value); }
};

template<typename S, typename V, typename M>
struct Parametrised2dTransformContext : public ParametrisedTransformContext<S, V, M, S*>
{
	using scalar_t = S;
	using vector_t = V;
	using matrix_t = M;
	using P = scalar_t*;

	void pushRotation(P parameter) { pushSingleParameterTransformation<RotateIJTransform<scalar_t, vector_t, 0, 1>>(parameter); }
	void pushRotation(scalar_t value) { pushStaticTransformation<RotateIJTransform<scalar_t, vector_t, 0, 1>>(value); }

	void pushRotationAboutAxis(P parameter, Vec3<scalar_t> const& axis) = delete;
	void pushRotationAboutAxis(scalar_t angle, Vec3<scalar_t> const& axis) = delete;
	void pushRotationAboutX(P parameter) = delete;
	void pushRotationAboutX(scalar_t angle) = delete;
	void pushRotationAboutY(P parameter) = delete;
	void pushRotationAboutY(scalar_t angle) = delete;
	void pushRotationAboutZ(P parameter) = delete;
	void pushRotationAboutZ(scalar_t angle) = delete;

	void pushScaleZ(P parameter) = delete;
	void pushScaleZ(scalar_t factor) = delete;

	void pushTranslationAlongZ(P parameter) = delete;
	void pushTranslationAlongZ(scalar_t value) = delete;
};