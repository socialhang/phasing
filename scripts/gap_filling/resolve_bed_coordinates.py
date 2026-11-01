#!/usr/bin/env python3
"""Convert integer arithmetic in a reviewed six-column fragment plan to BED.

Supports integer literals, addition, subtraction and parentheses.
"""
import ast
import sys


def integer(expression):
    node = ast.parse(expression, mode="eval").body

    def evaluate(item):
        if isinstance(item, ast.Constant) and type(item.value) is int:
            return item.value
        if isinstance(item, ast.BinOp) and isinstance(item.op, (ast.Add, ast.Sub)):
            left, right = evaluate(item.left), evaluate(item.right)
            return left + right if isinstance(item.op, ast.Add) else left - right
        if isinstance(item, ast.UnaryOp) and isinstance(item.op, (ast.UAdd, ast.USub)):
            value = evaluate(item.operand)
            return value if isinstance(item.op, ast.UAdd) else -value
        raise ValueError("Coordinates must use integers, +, - and parentheses")

    return evaluate(node)


with open(sys.argv[1], encoding="utf-8") as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        source, start, end, target, score, strand = line.split()
        start, end = integer(start), integer(end)
        if start < 0 or end <= start or strand not in ("+", "-"):
            raise ValueError("Invalid fragment interval: " + line.strip())
        print(source, start, end, target, score, strand, sep="\t")
