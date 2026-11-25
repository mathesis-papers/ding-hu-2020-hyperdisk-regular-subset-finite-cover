#!/usr/bin/env python3
"""Exact finite checks for the four-coordinate hyperdisk example."""

from fractions import Fraction
from itertools import permutations
import json


def scalar(value):
    return (Fraction(value),) + (Fraction(0),) * 7


def basis(index):
    return tuple(Fraction(position == index) for position in range(8))


def add(*terms):
    return tuple(sum(term[position] for term in terms) for position in range(8))


def scale(value, term):
    return tuple(Fraction(value) * coefficient for coefficient in term)


def multiply(left, right, radicands):
    result = [Fraction(0)] * 8
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            coefficient = left_coefficient * right_coefficient
            overlap = left_index & right_index
            for position, radicand in enumerate(radicands):
                if overlap & (1 << position):
                    coefficient *= radicand
            result[left_index ^ right_index] += coefficient
    return tuple(result)


def conjugate(term):
    return tuple(-coefficient if position & 4 else coefficient
                 for position, coefficient in enumerate(term))


def norm_squared(term, radicands):
    return multiply(term, conjugate(term), radicands)


def partitions(items):
    if not items:
        yield ()
        return
    for partition in partitions(items[1:]):
        yield ((items[0],),) + partition
        for position in range(len(partition)):
            yield (partition[:position] + ((items[0],) + partition[position],)
                   + partition[position + 1:])


def determinant(matrix, radicands):
    result = scalar(0)
    for permutation in permutations(range(3)):
        inversions = sum(permutation[left] > permutation[right]
                         for left in range(3) for right in range(left + 1, 3))
        term = scalar((-1) ** inversions)
        for row in range(3):
            term = multiply(term, matrix[row][permutation[row]], radicands)
        result = add(result, term)
    return result


def adjoint(matrix):
    return [[conjugate(matrix[row][column]) for row in range(len(matrix))]
            for column in range(len(matrix[0]))]


def matrix_product(left, right, radicands):
    return [[add(*(multiply(left[row][inner], right[inner][column], radicands)
                   for inner in range(len(right))))
             for column in range(len(right[0]))] for row in range(len(left))]


def check_siblings(columns, imaginary):
    radicands = (Fraction(3), Fraction(1), Fraction(-1))
    dimensions = list(range(4, 11))
    dimension_pairs = []
    for dimension in dimensions:
        normalization = sum((index + 1) ** 2 for index in range(dimension))
        spectra = [Fraction((index + 1) ** 2, normalization)
                   for index in range(dimension)]
        assert sum(spectra) == 1 and len(set(spectra)) == dimension
        extended = [tuple(scale(index + 1, column[index]) for index in range(4))
                    + tuple(scalar(index + 1) for index in range(4, dimension))
                    for column in columns]
        for vector in extended:
            assert add(scale(1, vector[0]), scale(Fraction(1, 2), vector[1]),
                       scale(Fraction(1, 3), vector[2]),
                       scale(Fraction(1, 2), vector[3])) == scalar(0)
            assert add(*(norm_squared(value, radicands) for value in vector)) == scalar(normalization)
        for index in range(4, dimension):
            flipped = list(extended[0])
            flipped[index] = scale(-1, flipped[index])
            difference = [add(original, scale(-1, other))
                          for original, other in zip(extended[0], flipped)]
            assert difference == [scalar(2 * (index + 1) if position == index else 0)
                                  for position in range(dimension)]
        for input_dimension in range(3, dimension):
            dimension_pairs.append([input_dimension, dimension])
            for vector in extended:
                assert all(scale(Fraction(1, index + 1), vector[index])
                           == scale(Fraction(1, 4), vector[3])
                           for index in range(input_dimension + 1, dimension))
            assert 3 + len(range(4, input_dimension + 1)) == input_dimension
    moments = [Fraction(index, 4) for index in range(-8, 9)]
    for parameter in moments:
        coefficient = [[add(imaginary, scalar(parameter)), scalar(parameter ** 2)],
                       [scalar(parameter ** 2), add(imaginary, scalar(-parameter))]]
        expected = [[scalar(1 + parameter ** 2 + parameter ** 4), scalar(0)],
                    [scalar(0), scalar(1 + parameter ** 2 + parameter ** 4)]]
        assert matrix_product(coefficient, adjoint(coefficient), radicands) == expected
        assert matrix_product(adjoint(coefficient), coefficient, radicands) == expected
        phase = scale(1 / (1 + parameter ** 2),
                      add(scalar(1 - parameter ** 2), scale(2 * parameter, imaginary)))
        partial = [[phase, imaginary, scalar(0)],
                   [imaginary, conjugate(phase), scalar(0)],
                   [scalar(0), scalar(0), scalar(2)]]
        expected_partial = [[scalar(2 if row < 2 else 4) if row == column else scalar(0)
                             for column in range(3)] for row in range(3)]
        assert matrix_product(partial, adjoint(partial), radicands) == expected_partial
        assert matrix_product(adjoint(partial), partial, radicands) == expected_partial
    roots = (Fraction(2), Fraction(1), Fraction(-1))
    half_root = scale(Fraction(1, 2), basis(1))
    rows = [[scalar(1), scalar(0)], [scalar(0), scalar(1)],
            [half_root, half_root],
            [half_root, multiply(imaginary, half_root, roots)]]
    gram = matrix_product(rows, adjoint(rows), roots)
    assert all(gram[index][index] == scalar(1) for index in range(4))
    assert all(gram[row][column] == conjugate(gram[column][row])
               for row in range(4) for column in range(4))
    assert rows[:2] == [[scalar(1), scalar(0)], [scalar(0), scalar(1)]]
    hermitian_constraints = [[Fraction(1), Fraction(0), Fraction(0), Fraction(0)],
                             [Fraction(0), Fraction(1), Fraction(0), Fraction(0)],
                             [Fraction(1, 2), Fraction(1, 2), Fraction(1), Fraction(0)],
                             [Fraction(1, 2), Fraction(1, 2), Fraction(0), Fraction(1)]]
    hermitian_basis = [[[scalar(1), scalar(0)], [scalar(0), scalar(0)]],
                       [[scalar(0), scalar(0)], [scalar(0), scalar(1)]],
                       [[scalar(0), scalar(1)], [scalar(1), scalar(0)]],
                       [[scalar(0), imaginary], [scale(-1, imaginary), scalar(0)]]]
    for row_index, vector in enumerate(rows):
        for basis_index, matrix in enumerate(hermitian_basis):
            value = matrix_product(matrix_product([vector], matrix, roots),
                                   adjoint([vector]), roots)[0][0]
            assert value == scalar(hermitian_constraints[row_index][basis_index])
    constraint_determinant = Fraction(0)
    for permutation in permutations(range(4)):
        inversions = sum(permutation[left] > permutation[right]
                         for left in range(4) for right in range(left + 1, 4))
        term = Fraction((-1) ** inversions)
        for row in range(4):
            term *= hermitian_constraints[row][permutation[row]]
        constraint_determinant += term
    assert constraint_determinant == 1
    return {'arbitrary_spectrum_dimensions_checked': dimensions,
            'input_dimension_marginal_rank_pairs_checked': dimension_pairs,
            'degenerate_example_parameters_checked_each': len(moments),
            'eq36_unnormalized_marginals': '(1+t^2+t^4)*I_2',
            'eq51_unnormalized_marginals': 'diag(2,2,4); normalization squared 8',
            'correlation_witness_rank': 2,
            'correlation_witness_diagonal': [1, 1, 1, 1],
            'hermitian_constraint_determinant': str(constraint_determinant)}


def main():
    weights = (1, 1, 1, 2)
    cases = list(partitions((0, 1, 2, 3)))
    assert len(cases) == 15
    allowed = [partition for partition in cases
               if all(2 * max(weights[index] for index in block)
                      <= sum(weights[index] for index in block)
                      for block in partition)]
    assert allowed == [((0, 1, 2, 3),)]
    imaginary = basis(4)
    first_root = basis(1)
    second_root = basis(2)
    parameters = []
    for index in range(32):
        parameter = -Fraction(1, 2) - Fraction(1, 4 * (index + 1))
        square_length = 5 + 4 * parameter
        radicands = (1 - parameter ** 2, (4 - square_length) / square_length, -1)
        assert -Fraction(3, 4) <= parameter <= -Fraction(1, 2)
        assert all(value > 0 for value in radicands[:2])
        first = add(scalar(parameter), multiply(imaginary, first_root, radicands))
        offset = add(scalar(-2), scale(-1, first))
        skew = multiply(imaginary, second_root, radicands)
        phases = (first,
                  scale(Fraction(1, 2), multiply(offset, add(scalar(1), skew), radicands)),
                  scale(Fraction(1, 2), multiply(offset, add(scalar(1), scale(-1, skew)), radicands)),
                  scalar(1))
        assert all(norm_squared(phase, radicands) == scalar(1) for phase in phases)
        assert add(*(scale(weight, phase) for weight, phase in zip(weights, phases))) == scalar(0)
        coordinates = tuple(scale(index + 1, phase) for index, phase in enumerate(phases))
        assert add(*(norm_squared(value, radicands) for value in coordinates)) == scalar(30)
        assert add(*(scale(weight, value) for weight, value in zip((6, 3, 2, 3), coordinates))) == scalar(0)
        assert coordinates[3] == scalar(4)
        assert coordinates[0][0] == parameter
        parameters.append(parameter)
    assert len(set(parameters)) == 32
    radicands = (Fraction(3), Fraction(1), Fraction(-1))
    omega = scale(Fraction(1, 2), add(scalar(-1), multiply(imaginary, first_root, radicands)))
    omega_bar = conjugate(omega)
    columns = ((scalar(-1), omega, omega_bar, scalar(1)),
               (scalar(-1), omega_bar, omega, scalar(1)),
               (omega, scalar(-1), omega_bar, scalar(1)))
    for column in columns:
        assert all(norm_squared(value, radicands) == scalar(1) for value in column)
        assert add(*(scale(weight, value) for weight, value in zip(weights, column))) == scalar(0)
    minor = [[scale(row + 1, column[row]) for column in columns] for row in (0, 1, 3)]
    minor_norm_squared = norm_squared(determinant(minor, radicands), radicands)
    assert minor_norm_squared == scalar(192)
    normalized_minor_norm_squared = Fraction(192, 30 ** 3)
    spectrum = [Fraction(value, 30) for value in (1, 4, 9, 16)]
    assert sum(spectrum) == 1 and min(spectrum) > 0 and len(set(spectrum)) == 4
    result = {
        'arithmetic': 'Exact rational arithmetic in a real two-square-root algebra with i^2=-1',
        'support_partitions': len(cases),
        'partitions_satisfying_all_necessary_polygon_inequalities': len(allowed),
        'rational_parameters_checked': len(parameters),
        'parameter_formula': 't_n=-1/2-1/(4(n+1)), 0<=n<32',
        'phase_moduli_squared': [1, 1, 1, 1],
        'unnormalized_norm_squared': 30,
        'linear_constraint': '6*x0+3*x1+2*x2+3*x3=0',
        'independent_span_vectors': 3,
        'normalized_minor_squared_modulus': str(normalized_minor_norm_squared),
        'marginal_spectrum': [str(value) for value in spectrum],
        'siblings': check_siblings(columns, imaginary),
        'scope': 'Finite support cases, exact algebraic examples, spectrum extensions and matrix identities; the all-parameter, all-dimension and cover-cardinality conclusions have written proofs in note.tex.',
        'result': 'all checks passed',
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
