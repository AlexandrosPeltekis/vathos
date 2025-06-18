#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include <functional>
#include <unordered_set>

// Define NodePtr early
struct Node;
using NodePtr = std::shared_ptr<Node>;

struct Node {
    double val;
    double grad = 0.0;
    std::function<void()> backward;
    std::vector<NodePtr> inputs;
};

struct rVar {
    NodePtr node;

    // Constructor
    rVar(double val) : node(std::make_shared<Node>()) {
        node->val = val;
    }

    // Value and gradient access
    double val() const { return node->val; }
    double grad() const { return node->grad; }

    // Backpropagation
    void backward() {
        node->grad = 1.0;
        std::vector<NodePtr> stack = { node };
        std::unordered_set<NodePtr> visited;

        while (!stack.empty()) {
            NodePtr n = stack.back();
            stack.pop_back();

            if (visited.count(n)) continue;
            visited.insert(n);

            if (n->backward) n->backward();

            for (const NodePtr& input : n->inputs) {
                if (input) stack.push_back(input);
            }
        }
    }

    // Operator Overloads
    friend rVar operator+(const rVar& a, const rVar& b) {
        rVar out(a.val() + b.val());
        out.node->inputs = { a.node, b.node };
        out.node->backward = [a, b, self = out.node]() {
            std::cout << "Running backward for +\n";
            a.node->grad += self->grad;
            b.node->grad += self->grad;
            };
        return out;
    }

    friend rVar operator*(const rVar& a, const rVar& b) {
        rVar out(a.val() * b.val());
        out.node->inputs = { a.node, b.node };
        out.node->backward = [a, b, self = out.node]() {
            a.node->grad += b.val() * self->grad;
            b.node->grad += a.val() * self->grad;
            };
        return out;
    }

    friend rVar operator-(const rVar& a, const rVar& b) {
        rVar out(a.val() - b.val());
        out.node->inputs = { a.node, b.node };
        out.node->backward = [a, b, self = out.node]() {
            a.node->grad += self->grad;
            b.node->grad -= self->grad;
            };
        return out;
    }

    friend rVar operator/(const rVar& a, const rVar& b) {
        rVar out(a.val() / b.val());
        out.node->inputs = { a.node, b.node };
        out.node->backward = [a, b, self = out.node]() {
            a.node->grad += (1.0 / b.val()) * self->grad;
            b.node->grad -= (a.val() / (b.val() * b.val())) * self->grad;
            };
        return out;
    }

    friend rVar sin(const rVar& a) {
        rVar out(std::sin(a.val()));
        out.node->inputs = { a.node };
        out.node->backward = [a, self = out.node]() {
            a.node->grad += std::cos(a.val()) * self->grad;
            };
        return out;
    }

    friend rVar cos(const rVar& a) {
        rVar out(std::cos(a.val()));
        out.node->inputs = { a.node };
        out.node->backward = [a, self = out.node]() {
            a.node->grad -= std::sin(a.val()) * self->grad;
            };
        return out;
    }

    friend rVar exp(const rVar& a) {
        rVar out(std::exp(a.val()));
        out.node->inputs = { a.node };
        out.node->backward = [a, self = out.node]() {
            a.node->grad += self->val * self->grad;
            };
        return out;
    }

    friend rVar pow(const rVar& a, double p) {
        rVar out(std::pow(a.val(), p));
        out.node->inputs = { a.node };
        out.node->backward = [a, p, self = out.node]() {
            a.node->grad += p * std::pow(a.val(), p - 1) * self->grad;
            };
        return out;
    }
};
