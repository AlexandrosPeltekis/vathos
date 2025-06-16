#pragma once
#include <vector>
#include <memory>
#include <functional>
#include <cassert>
#include <array>
#include <cmath>
#include <string>
#include <iostream>
#include <unordered_set>

// Reverse-mode autodiff variable
class RevVar {
public:
    double val = 0.0;   // Value of the variable
    double grad = 0.0;  // Gradient (to be computed in backward pass)

    // Graph info
    std::string op = "input"; // Operation name
    std::vector<std::weak_ptr<RevVar>> parents; // Parent nodes
    int id = -1; // Unique node id

    // Each node in the computation graph stores how to propagate gradients to its parents
    using BackwardFunc = std::function<void()>;
    std::vector<BackwardFunc> backward_ops;

    // For memory management, keep all variables alive
    static std::vector<std::shared_ptr<RevVar>>& tape() {
        static std::vector<std::shared_ptr<RevVar>> t;
        return t;
    }

    // For unique node IDs
    static int& next_id() {
        static int id = 0;
        return id;
    }

    // Constructors
    RevVar() : val(0.0), id(next_id()++) { tape().emplace_back(this); }
    explicit RevVar(double v) : val(v), id(next_id()++) { tape().emplace_back(this); }

    // Factory for shared_ptr
    static std::shared_ptr<RevVar> create(double v, const std::string& op = "input", const std::vector<std::shared_ptr<RevVar>>& parents = {}) {
        auto ptr = std::make_shared<RevVar>(v);
        ptr->op = op;
        for (const auto& p : parents) {
            ptr->parents.push_back(p);
        }
        tape().push_back(ptr);
        return ptr;
    }

    // Operator overloads
    friend std::shared_ptr<RevVar> operator+(const std::shared_ptr<RevVar>& a, const std::shared_ptr<RevVar>& b) {
        auto out = RevVar::create(a->val + b->val, "+", {a, b});
        out->backward_ops.push_back([a, out]() { a->grad += out->grad; });
        out->backward_ops.push_back([b, out]() { b->grad += out->grad; });
        return out;
    }
    friend std::shared_ptr<RevVar> operator-(const std::shared_ptr<RevVar>& a, const std::shared_ptr<RevVar>& b) {
        auto out = RevVar::create(a->val - b->val, "-", {a, b});
        out->backward_ops.push_back([a, out]() { a->grad += out->grad; });
        out->backward_ops.push_back([b, out]() { b->grad -= out->grad; });
        return out;
    }
    friend std::shared_ptr<RevVar> operator*(const std::shared_ptr<RevVar>& a, const std::shared_ptr<RevVar>& b) {
        auto out = RevVar::create(a->val * b->val, "*", {a, b});
        out->backward_ops.push_back([a, b, out]() { a->grad += b->val * out->grad; });
        out->backward_ops.push_back([a, b, out]() { b->grad += a->val * out->grad; });
        return out;
    }
    friend std::shared_ptr<RevVar> operator/(const std::shared_ptr<RevVar>& a, const std::shared_ptr<RevVar>& b) {
        auto out = RevVar::create(a->val / b->val, "/", {a, b});
        out->backward_ops.push_back([a, b, out]() { a->grad += out->grad / b->val; });
        out->backward_ops.push_back([a, b, out]() { b->grad -= (a->val / (b->val * b->val)) * out->grad; });
        return out;
    }

    // Overloads for scalar (double) on right
    friend std::shared_ptr<RevVar> operator+(const std::shared_ptr<RevVar>& a, double b) {
        return RevVar::create(a->val + b, "+", {a});
    }
    friend std::shared_ptr<RevVar> operator-(const std::shared_ptr<RevVar>& a, double b) {
        return RevVar::create(a->val - b, "-", {a});
    }
    friend std::shared_ptr<RevVar> operator*(const std::shared_ptr<RevVar>& a, double b) {
        auto out = RevVar::create(a->val * b, "*", {a});
        out->backward_ops.push_back([a, b, out]() { a->grad += b * out->grad; });
        return out;
    }
    friend std::shared_ptr<RevVar> operator/(const std::shared_ptr<RevVar>& a, double b) {
        auto out = RevVar::create(a->val / b, "/", {a});
        out->backward_ops.push_back([a, b, out]() { a->grad += out->grad / b; });
        return out;
    }

    // Overloads for scalar (double) on left
    friend std::shared_ptr<RevVar> operator+(double a, const std::shared_ptr<RevVar>& b) {
        return RevVar::create(a + b->val, "+", {b});
    }
    friend std::shared_ptr<RevVar> operator-(double a, const std::shared_ptr<RevVar>& b) {
        auto out = RevVar::create(a - b->val, "-", {b});
        out->backward_ops.push_back([b, out]() { b->grad -= out->grad; });
        return out;
    }
    friend std::shared_ptr<RevVar> operator*(double a, const std::shared_ptr<RevVar>& b) {
        auto out = RevVar::create(a * b->val, "*", {b});
        out->backward_ops.push_back([b, a, out]() { b->grad += a * out->grad; });
        return out;
    }
    friend std::shared_ptr<RevVar> operator/(double a, const std::shared_ptr<RevVar>& b) {
        auto out = RevVar::create(a / b->val, "/", {b});
        out->backward_ops.push_back([b, a, out]() { b->grad -= (a / (b->val * b->val)) * out->grad; });
        return out;
    }

    // Elementary functions
    friend std::shared_ptr<RevVar> sin(const std::shared_ptr<RevVar>& a) {
        auto out = RevVar::create(std::sin(a->val), "sin", {a});
        out->backward_ops.push_back([a, out]() { a->grad += std::cos(a->val) * out->grad; });
        return out;
    }
    friend std::shared_ptr<RevVar> cos(const std::shared_ptr<RevVar>& a) {
        auto out = RevVar::create(std::cos(a->val), "cos", {a});
        out->backward_ops.push_back([a, out]() { a->grad -= std::sin(a->val) * out->grad; });
        return out;
    }
    friend std::shared_ptr<RevVar> exp(const std::shared_ptr<RevVar>& a) {
        auto out = RevVar::create(std::exp(a->val), "exp", {a});
        out->backward_ops.push_back([a, out]() { a->grad += out->val * out->grad; });
        return out;
    }
    friend std::shared_ptr<RevVar> log(const std::shared_ptr<RevVar>& a) {
        auto out = RevVar::create(std::log(a->val), "log", {a});
        out->backward_ops.push_back([a, out]() { a->grad += (1.0 / a->val) * out->grad; });
        return out;
    }

    // Backward pass: call on output variable
    void backward() {
        grad = 1.0;
        for (auto it = tape().rbegin(); it != tape().rend(); ++it) {
            for (auto& op : (*it)->backward_ops) {
                op();
            }
        }
    }

    // Reset tape and gradients
    static void clear_tape() {
        for (auto& v : tape()) v->grad = 0.0;
        tape().clear();
        next_id() = 0;
    }

    // Print the computation graph (DFS)
    void print_graph(std::ostream& os = std::cout) const {
        std::unordered_set<int> visited;
        print_graph_impl(os, visited, 0);
    }

private:
    void print_graph_impl(std::ostream& os, std::unordered_set<int>& visited, int depth) const {
        if (visited.count(id)) return;
        visited.insert(id);
        for (int i = 0; i < depth; ++i) os << "  ";
        os << "Node " << id << ": op=" << op << ", val=" << val << ", grad=" << grad << "\n";
        for (const auto& wp : parents) {
            if (auto p = wp.lock()) {
                p->print_graph_impl(os, visited, depth + 1);
            }
        }
    }
};

// Reverse-mode autodiff variable for 3 variables (x, y, z) with Hessian
class RevVar3 {
public:
    double val = 0.0; // Value of the variable

    // Gradient: df/dx, df/dy, df/dz
    std::array<double, 3> grad{0.0, 0.0, 0.0};

    // Hessian: d²f/dx², d²f/dxdy, d²f/dxdz, d²f/dy², d²f/dydz, d²f/dz²
    // Stored as: [xx, xy, xz, yy, yz, zz]
    std::array<double, 6> hess{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    using BackwardFunc = std::function<void()>;
    std::vector<BackwardFunc> backward_ops;

    static std::vector<std::shared_ptr<RevVar3>>& tape() {
        static std::vector<std::shared_ptr<RevVar3>> t;
        return t;
    }

    RevVar3() : val(0.0) { tape().emplace_back(this); }
    explicit RevVar3(double v) : val(v) { tape().emplace_back(this); }

    static std::shared_ptr<RevVar3> create(double v) {
        auto ptr = std::make_shared<RevVar3>(v);
        tape().push_back(ptr);
        return ptr;
    }

    // Helper: index for Hessian
    static constexpr int idx(int i, int j) {
        // (0,0)->0, (0,1)->1, (0,2)->2, (1,1)->3, (1,2)->4, (2,2)->5
        if (i == 0 && j == 0) return 0;
        if ((i == 0 && j == 1) || (i == 1 && j == 0)) return 1;
        if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 2;
        if (i == 1 && j == 1) return 3;
        if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 4;
        if (i == 2 && j == 2) return 5;
        assert(false);
        return -1;
    }

    // Operator overloads
    friend std::shared_ptr<RevVar3> operator+(const std::shared_ptr<RevVar3>& a, const std::shared_ptr<RevVar3>& b) {
        auto out = RevVar3::create(a->val + b->val);
        out->backward_ops.push_back([a, out]() {
            for (int i = 0; i < 3; ++i) a->grad[i] += out->grad[i];
            for (int i = 0; i < 6; ++i) a->hess[i] += out->hess[i];
        });
        out->backward_ops.push_back([b, out]() {
            for (int i = 0; i < 3; ++i) b->grad[i] += out->grad[i];
            for (int i = 0; i < 6; ++i) b->hess[i] += out->hess[i];
        });
        return out;
    }
    friend std::shared_ptr<RevVar3> operator-(const std::shared_ptr<RevVar3>& a, const std::shared_ptr<RevVar3>& b) {
        auto out = RevVar3::create(a->val - b->val);
        out->backward_ops.push_back([a, out]() {
            for (int i = 0; i < 3; ++i) a->grad[i] += out->grad[i];
            for (int i = 0; i < 6; ++i) a->hess[i] += out->hess[i];
        });
        out->backward_ops.push_back([b, out]() {
            for (int i = 0; i < 3; ++i) b->grad[i] -= out->grad[i];
            for (int i = 0; i < 6; ++i) b->hess[i] -= out->hess[i];
        });
        return out;
    }
    friend std::shared_ptr<RevVar3> operator*(const std::shared_ptr<RevVar3>& a, const std::shared_ptr<RevVar3>& b) {
        auto out = RevVar3::create(a->val * b->val);
        out->backward_ops.push_back([a, b, out]() {
            // Gradient
            for (int i = 0; i < 3; ++i)
                a->grad[i] += b->val * out->grad[i];
            // Hessian
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j <= i; ++j)
                    a->hess[idx(i, j)] += b->val * out->hess[idx(i, j)];
            // Mixed terms
            for (int i = 0; i < 3; ++i)
                a->grad[i] += b->grad[i] * out->val;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j <= i; ++j)
                    a->hess[idx(i, j)] += b->grad[i] * out->grad[j] + b->grad[j] * out->grad[i];
        });
        out->backward_ops.push_back([a, b, out]() {
            for (int i = 0; i < 3; ++i)
                b->grad[i] += a->val * out->grad[i];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j <= i; ++j)
                    b->hess[idx(i, j)] += a->val * out->hess[idx(i, j)];
            for (int i = 0; i < 3; ++i)
                b->grad[i] += a->grad[i] * out->val;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j <= i; ++j)
                    b->hess[idx(i, j)] += a->grad[i] * out->grad[j] + a->grad[j] * out->grad[i];
        });
        return out;
    }
    // For brevity, division and scalar overloads are omitted but can be implemented similarly.

    // Overloads for scalar (double) on right
    friend std::shared_ptr<RevVar3> operator+(const std::shared_ptr<RevVar3>& a, double b) {
        return RevVar3::create(a->val + b);
    }
    friend std::shared_ptr<RevVar3> operator-(const std::shared_ptr<RevVar3>& a, double b) {
        return RevVar3::create(a->val - b);
    }
    friend std::shared_ptr<RevVar3> operator*(const std::shared_ptr<RevVar3>& a, double b) {
        auto out = RevVar3::create(a->val * b);
        out->backward_ops.push_back([a, b, out]() {
            for (int i = 0; i < 3; ++i) a->grad[i] += b * out->grad[i];
            for (int i = 0; i < 6; ++i) a->hess[i] += b * out->hess[i];
        });
        return out;
    }

    // Overloads for scalar (double) on left
    friend std::shared_ptr<RevVar3> operator+(double a, const std::shared_ptr<RevVar3>& b) {
        return RevVar3::create(a + b->val);
    }
    friend std::shared_ptr<RevVar3> operator-(double a, const std::shared_ptr<RevVar3>& b) {
        auto out = RevVar3::create(a - b->val);
        out->backward_ops.push_back([b, out]() {
            for (int i = 0; i < 3; ++i) b->grad[i] -= out->grad[i];
            for (int i = 0; i < 6; ++i) b->hess[i] -= out->hess[i];
        });
        return out;
    }
    friend std::shared_ptr<RevVar3> operator*(double a, const std::shared_ptr<RevVar3>& b) {
        auto out = RevVar3::create(a * b->val);
        out->backward_ops.push_back([b, a, out]() {
            for (int i = 0; i < 3; ++i) b->grad[i] += a * out->grad[i];
            for (int i = 0; i < 6; ++i) b->hess[i] += a * out->hess[i];
        });
        return out;
    }

    // Elementary functions (example: sin)
    friend std::shared_ptr<RevVar3> sin(const std::shared_ptr<RevVar3>& a) {
        auto out = RevVar3::create(std::sin(a->val));
        out->backward_ops.push_back([a, out]() {
            double c = std::cos(a->val);
            double s = std::sin(a->val);
            for (int i = 0; i < 3; ++i)
                a->grad[i] += c * out->grad[i];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j <= i; ++j)
                    a->hess[idx(i, j)] += -s * out->grad[i] * out->grad[j] + c * out->hess[idx(i, j)];
        });
        return out;
    }

    // Backward pass: call on output variable
    void backward(int wrt = 0) {
        grad[wrt] = 1.0;
        for (auto it = tape().rbegin(); it != tape().rend(); ++it) {
            for (auto& op : (*it)->backward_ops) {
                op();
            }
        }
    }

    // Reset tape and derivatives
    static void clear_tape() {
        for (auto& v : tape()) {
            v->grad = {0.0, 0.0, 0.0};
            v->hess = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        }
        tape().clear();
    }
};