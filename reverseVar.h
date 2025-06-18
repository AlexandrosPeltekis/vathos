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

	//Zero gradient
    void zero_grad() {
        std::vector<NodePtr> stack = { node };
        std::unordered_set<NodePtr> visited;
        while (!stack.empty()) {
            NodePtr n = stack.back();
            stack.pop_back();
            if (visited.count(n)) continue;
            visited.insert(n);
            n->grad = 0.0;
            for (const NodePtr& input : n->inputs) {
                if (input) stack.push_back(input);
            }
        }
    }

    // Operator Overloads
    friend rVar operator+(const rVar& a, const rVar& b) {
        rVar out(a.val() + b.val());
        out.node->inputs = { a.node, b.node };
        Node* a_ptr = a.node.get();
        Node* b_ptr = b.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, b_ptr, self]() {
            a_ptr->grad += self->grad;
            b_ptr->grad += self->grad;
            };
        return out;
    }

    friend rVar operator*(const rVar& a, const rVar& b) {
        rVar out(a.val() * b.val());
        out.node->inputs = { a.node, b.node };
        Node* a_ptr = a.node.get();
        Node* b_ptr = b.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, b_ptr, self]() {
            a_ptr->grad += b_ptr->val * self->grad;
            b_ptr->grad += a_ptr->val * self->grad;
            };
        return out;
    }

    friend rVar operator-(const rVar& a, const rVar& b) {
        rVar out(a.val() - b.val());
        out.node->inputs = { a.node, b.node };
        Node* a_ptr = a.node.get();
        Node* b_ptr = b.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, b_ptr, self]() {
            a_ptr->grad += self->grad;
            b_ptr->grad -= self->grad;
            };
        return out;
    }

    friend rVar operator/(const rVar& a, const rVar& b) {
        rVar out(a.val() / b.val());
        out.node->inputs = { a.node, b.node };
        Node* a_ptr = a.node.get();
        Node* b_ptr = b.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, b_ptr, self]() {
            a_ptr->grad += (1.0 / b_ptr->val) * self->grad;
            b_ptr->grad -= (a_ptr->val / (b_ptr->val * b_ptr->val)) * self->grad;
            };
        return out;
    }

    friend rVar sin(const rVar& a) {
        rVar out(std::sin(a.val()));
        out.node->inputs = { a.node };
        Node* a_ptr = a.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, self ]() {
            a_ptr->grad += std::cos(a_ptr->val) * self->grad;
            };
        return out;
    }

    friend rVar cos(const rVar& a) {
        rVar out(std::cos(a.val()));
        out.node->inputs = { a.node };
        Node* a_ptr = a.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, self ]() {
            a_ptr->grad -= std::sin(a_ptr->val) * self->grad;
            };
        return out;
    }

    friend rVar tan(const rVar& a) {
        rVar out(std::tan(a.val()));
        out.node->inputs = { a.node };
        Node* a_ptr = a.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, self ]() {
            a_ptr->grad += (1.0 / std::cos(a_ptr->val) / std::cos(a_ptr->val)) * self->grad;
            };
        return out;
	}

    friend rVar exp(const rVar& a) {
        rVar out(std::exp(a.val()));
        out.node->inputs = { a.node };
        Node* a_ptr = a.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, self ]() {
            a_ptr->grad += self->val * self->grad;
            };
        return out;
    }

    friend rVar pow(const rVar& a, double p) {
        rVar out(std::pow(a.val(), p));
        out.node->inputs = { a.node };
        Node* a_ptr = a.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, p, self]() {
            a_ptr->grad += p * std::pow(a_ptr->val, p - 1) * self->grad;
            };
        return out;
    }

    friend rVar pow(double base, const rVar& a) {
        rVar out(std::pow(base, a.val()));
        out.node->inputs = { a.node };
        Node* a_ptr = a.node.get();
        Node* self = out.node.get();
        out.node->backward = [base, a_ptr, self]() {
            a_ptr->grad += std::log(base) * std::pow(base, a_ptr->val) * self->grad;
            };
        return out;
	}

    friend rVar pow(const rVar& a, const rVar& b) {
        rVar out(std::pow(a.val(), b.val()));
        Node* a_ptr = a.node.get();
        Node* b_ptr = b.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, b_ptr, self]() {
            a_ptr->grad += b_ptr->val * std::pow(a_ptr->val, b_ptr->val - 1) * self->grad;
            b_ptr->grad += std::log(a_ptr->val) * std::pow(a_ptr->val, b_ptr->val) * self->grad;
            };
        return out;
	}

    friend rVar sqrt(const rVar& a) {
        rVar out(std::sqrt(a.val()));
        out.node->inputs = { a.node };
        Node* a_ptr = a.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, self ]() {
            a_ptr->grad += (0.5 / std::sqrt(a_ptr->val)) * self->grad;
            };
        return out;
	}

    friend rVar log(const rVar& a) {
        rVar out(std::log(a.val()));
        out.node->inputs = { a.node };
        Node* a_ptr = a.node.get();
        Node* self = out.node.get();
        out.node->backward = [a_ptr, self ]() {
            a_ptr->grad += (1.0 / a_ptr->val) * self->grad;
            };
        return out;
	}


};
