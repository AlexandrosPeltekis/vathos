#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include <functional>
#include <unordered_set>
#include <complex>
#include <type_traits>

/*
    rVar<T> - Mathematical Guide

    This class represents a variable for reverse-mode automatic differentiation (AD).
    It tracks the value of a function and enables efficient computation of gradients
    (first derivatives) with respect to its inputs via backpropagation.

    Members:
      val   = f(x)        // function value
      grad  = df/dx       // gradient (computed after calling backward())

    Constructors:
      rVar()
        - Initializes value to zero.
        - Represents the zero function.

      rVar(val)
        - Initializes value to 'val', gradient zero.
        - Represents a constant or variable.

    Usage for Automatic Differentiation:
      - To represent a variable x: rVar<T> x(x0);
      - To represent a constant:  rVar<T> c(c0);

    Backpropagation:
      - After constructing a computation graph using rVar operations,
        call f.backward() on the output variable f to compute gradients.
      - Each input variable's grad() will then contain df/dx for that variable.

    Example:
      rVar<double> x(1.0), y(2.0);
      rVar<double> f = x * x * y + sin(x);
      f.backward();
      x.grad() // now holds df/dx
      y.grad() // holds df/dy

    Notes:
      - The computation graph is built dynamically as you perform operations.
      - Use zero_grad() to reset all gradients in the graph before a new backward pass.
      - Supports scalar and vector reverse-mode AD, and can be composed with dual numbers for higher-order derivatives.

    Supported Operations:
      - Arithmetic: +, -, *, /
      - Math functions: sin, cos, tan, exp, log, sqrt, pow
      - All operations propagate the computation graph for correct gradient calculation.
*/

// Node definition
template<typename T>
struct Node;
template<typename T>
using NodePtr = std::shared_ptr<Node<T>>;

template<typename T>
struct Node {
    T val;
    T grad = T(0);
    std::function<void()> backward;
    std::vector<NodePtr<T>> inputs;
};

// rVar: reverse-mode autodiff variable
template<typename T>
struct rVar {
    using value_type = T;
    NodePtr<T> node;

    // Constructors
    rVar() : node(std::make_shared<Node<T>>()) { node->val = T(0.0); }
    rVar(const T& val) : node(std::make_shared<Node<T>>()) { 
        static_assert(!std::is_integral_v<T>, "rVar: Input (T) must not be an integral type. Use float, double, std::complex float or double, Dual1 or Dual2.");
        node->val = static_cast<T>(val); 
    }

    // Conversion constructor for type flexibility (still allowed, but doesn't connect graphs)
    template<typename U, typename = std::enable_if_t<!std::is_same_v<T, U>>>
    explicit rVar(const rVar<U>& other) : node(std::make_shared<Node<T>>()) {
        node->val = static_cast<T>(other.val());
    }

    T val() const { return node->val; }
    T grad() const { return node->grad; }

    // Backpropagation
    void backward() {
        node->grad = T(1);
        std::vector<NodePtr<T>> stack = { node };
        std::unordered_set<NodePtr<T>> visited;
        while (!stack.empty()) {
            NodePtr<T> n = stack.back();
            stack.pop_back();
            if (visited.count(n)) continue;
            visited.insert(n);
            if (n->backward) n->backward();
            for (const NodePtr<T>& input : n->inputs) {
                if (input) stack.push_back(input);
            }
        }
    }

    void zero_grad() {
        std::vector<NodePtr<T>> stack = { node };
        std::unordered_set<NodePtr<T>> visited;
        while (!stack.empty()) {
            NodePtr<T> n = stack.back();
            stack.pop_back();
            if (visited.count(n)) continue;
            visited.insert(n);
            n->grad = T(0);
            for (const NodePtr<T>& input : n->inputs) {
                if (input) stack.push_back(input);
            }
        }
    }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, const rVar& x) {
        os << x.val() << " (grad=" << x.grad() << ")";
        return os;
    }

    // --- rVar <op> rVar ---
    template<typename U>
    auto operator+(const rVar<U>& o) const {
        static_assert(std::is_same_v<T, U>, "rVar: Both operands must have the same type. Upcast your variables before building the graph.");
        rVar out(this->val() + o.val());
        out.node->inputs = { this->node, o.node };
        Node<T>* a_ptr = this->node.get();
        Node<T>* b_ptr = o.node.get();
        Node<T>* self = out.node.get();
        out.node->backward = [a_ptr, b_ptr, self]() {
            a_ptr->grad += self->grad;
            b_ptr->grad += self->grad;
            };
        return out;
    }

    template<typename U>
    auto operator-(const rVar<U>& o) const {
        static_assert(std::is_same_v<T, U>, "rVar: Both operands must have the same type. Upcast your variables before building the graph.");
        rVar out(this->val() - o.val());
        out.node->inputs = { this->node, o.node };
        Node<T>* a_ptr = this->node.get();
        Node<T>* b_ptr = o.node.get();
        Node<T>* self = out.node.get();
        out.node->backward = [a_ptr, b_ptr, self]() {
            a_ptr->grad += self->grad;
            b_ptr->grad -= self->grad;
            };
        return out;
    }

    template<typename U>
    auto operator*(const rVar<U>& o) const {
        static_assert(std::is_same_v<T, U>, "rVar: Both operands must have the same type. Upcast your variables before building the graph.");
        rVar out(this->val() * o.val());
        out.node->inputs = { this->node, o.node };
        Node<T>* a_ptr = this->node.get();
        Node<T>* b_ptr = o.node.get();
        Node<T>* self = out.node.get();
        out.node->backward = [a_ptr, b_ptr, self]() {
            a_ptr->grad += b_ptr->val * self->grad;
            b_ptr->grad += a_ptr->val * self->grad;
            };
        return out;
    }

    template<typename U>
    auto operator/(const rVar<U>& o) const {
        static_assert(std::is_same_v<T, U>, "rVar: Both operands must have the same type. Upcast your variables before building the graph.");
        rVar out(this->val() / o.val());
        out.node->inputs = { this->node, o.node };
        Node<T>* a_ptr = this->node.get();
        Node<T>* b_ptr = o.node.get();
        Node<T>* self = out.node.get();
        out.node->backward = [a_ptr, b_ptr, self]() {
            a_ptr->grad += T(1) / b_ptr->val * self->grad;
            b_ptr->grad -= (a_ptr->val / (b_ptr->val * b_ptr->val)) * self->grad;
            };
        return out;
    }

    // --- rVar <op> scalar ---
    template<typename U>
    auto operator+(U o) const {
        return *this + rVar<T>(static_cast<T>(o));
    }
    template<typename U>
    auto operator-(U o) const {
        return *this - rVar<T>(static_cast<T>(o));
    }
    template<typename U>
    auto operator*(U o) const {
        return *this * rVar<T>(static_cast<T>(o));
    }
    template<typename U>
    auto operator/(U o) const {
        return *this / rVar<T>(static_cast<T>(o));
    }
};

// Deduction guide for rVar
template<typename T>
rVar(T) -> rVar<T>;

// --- Scalar <op> rVar ---
template<typename T, typename U>
auto operator+(U lhs, const rVar<T>& rhs) { return rVar<T>(static_cast<T>(lhs)) + rhs; }
template<typename T, typename U>
auto operator-(U lhs, const rVar<T>& rhs) { return rVar<T>(static_cast<T>(lhs)) - rhs; }
template<typename T, typename U>
auto operator*(U lhs, const rVar<T>& rhs) { return rVar<T>(static_cast<T>(lhs)) * rhs; }
template<typename T, typename U>
auto operator/(U lhs, const rVar<T>& rhs) { return rVar<T>(static_cast<T>(lhs)) / rhs; }

// --- Math functions (ADL-friendly) ---
template<typename T>
auto sin(const rVar<T>& a) {
    using std::sin;
    rVar<T> out(sin(a.val()));
    out.node->inputs = { a.node };
    Node<T>* a_ptr = a.node.get();
    Node<T>* self = out.node.get();
    out.node->backward = [a_ptr, self]() {
        using std::cos;
        a_ptr->grad += cos(a_ptr->val) * self->grad;
        };
    return out;
}
template<typename T>
auto cos(const rVar<T>& a) {
    using std::cos;
    rVar<T> out(cos(a.val()));
    out.node->inputs = { a.node };
    Node<T>* a_ptr = a.node.get();
    Node<T>* self = out.node.get();
    out.node->backward = [a_ptr, self]() {
        using std::sin;
        a_ptr->grad -= sin(a_ptr->val) * self->grad;
        };
    return out;
}
template<typename T>
auto tan(const rVar<T>& a) {
    using std::tan;
    using std::cos;
    rVar<T> out(tan(a.val()));
    out.node->inputs = { a.node };
    Node<T>* a_ptr = a.node.get();
    Node<T>* self = out.node.get();
    out.node->backward = [a_ptr, self]() {
        a_ptr->grad += (T(1) / (cos(a_ptr->val) * cos(a_ptr->val))) * self->grad;
        };
    return out;
}
template<typename T>
auto exp(const rVar<T>& a) {
    using std::exp;
    rVar<T> out(exp(a.val()));
    out.node->inputs = { a.node };
    Node<T>* a_ptr = a.node.get();
    Node<T>* self = out.node.get();
    out.node->backward = [a_ptr, self]() {
        a_ptr->grad += self->val * self->grad;
        };
    return out;
}
template<typename T>
auto log(const rVar<T>& a) {
    using std::log;
    rVar<T> out(log(a.val()));
    out.node->inputs = { a.node };
    Node<T>* a_ptr = a.node.get();
    Node<T>* self = out.node.get();
    out.node->backward = [a_ptr, self]() {
        a_ptr->grad += (T(1) / a_ptr->val) * self->grad;
        };
    return out;
}
template<typename T>
auto sqrt(const rVar<T>& a) {
    using std::sqrt;
    rVar<T> out(sqrt(a.val()));
    out.node->inputs = { a.node };
    Node<T>* a_ptr = a.node.get();
    Node<T>* self = out.node.get();
    out.node->backward = [a_ptr, self]() {
        a_ptr->grad += (T(0.5) / sqrt(a_ptr->val)) * self->grad;
        };
    return out;
}
template<typename T>
auto pow(const rVar<T>& a, T p) {
    using std::pow;
    rVar<T> out(pow(a.val(), p));
    out.node->inputs = { a.node };
    Node<T>* a_ptr = a.node.get();
    Node<T>* self = out.node.get();
    out.node->backward = [a_ptr, p, self]() {
        using std::pow;
        a_ptr->grad += p * pow(a_ptr->val, p - T(1)) * self->grad;
        };
    return out;
}
template<typename T>
auto pow(T base, const rVar<T>& a) {
    using std::pow;
    using std::log;
    rVar<T> out(pow(base, a.val()));
    out.node->inputs = { a.node };
    Node<T>* a_ptr = a.node.get();
    Node<T>* self = out.node.get();
    out.node->backward = [base, a_ptr, self]() {
        using std::pow;
        using std::log;
        a_ptr->grad += log(base) * pow(base, a_ptr->val) * self->grad;
        };
    return out;
}
template<typename T>
auto pow(const rVar<T>& a, const rVar<T>& b) {
    using std::pow;
    using std::log;
    rVar<T> out(pow(a.val(), b.val()));
    out.node->inputs = { a.node, b.node };
    Node<T>* a_ptr = a.node.get();
    Node<T>* b_ptr = b.node.get();
    Node<T>* self = out.node.get();
    out.node->backward = [a_ptr, b_ptr, self]() {
        using std::pow;
        using std::log;
        a_ptr->grad += b_ptr->val * pow(a_ptr->val, b_ptr->val - T(1)) * self->grad;
        b_ptr->grad += log(a_ptr->val) * pow(a_ptr->val, b_ptr->val) * self->grad;
        };
    return out;
}