"""Two-qubit ZY-rotation gate."""

import math
import numpy as np
from typing import Optional
from qiskit.circuit.gate import Gate
from qiskit.circuit.quantumregister import QuantumRegister
from qiskit.circuit.parameterexpression import ParameterValueType


class RZYGate(Gate):
    r"""A parametric 2-qubit :math:`Z \otimes Y` interaction (rotation about ZY).

    **Circuit Symbol:**

    .. parsed-literal::

             ┌─────────┐
        q_0: ┤0        ├
             │  Rzy(θ) │
        q_1: ┤1        ├
             └─────────┘

    **Matrix Representation:**

    .. math::

        \newcommand{\rotationangle}{\frac{\theta}{2}}

        R_{ZY}(\theta)\ q_0, q_1 = \exp\left(-i \frac{\theta}{2} Y{\otimes}Z\right) =
            \begin{pmatrix}
                \cos\left(\rotationangle\right) & 0 & -\sin\left(\rotationangle\right) & 0 \\
                0 & \cos\left(\rotationangle\right) & 0 & \sin\left(\rotationangle\right) \\
                \sin\left(\rotationangle\right) & 0 & \cos\left(\rotationangle\right) & 0 \\
                0 & -\sin\left(\rotationangle\right) & 0 & \cos\left(\rotationangle\right)
            \end{pmatrix}

    .. note::

        In Qiskit's convention, higher qubit indices are more significant
        (little endian convention). In the above example we apply the gate
        on (q_0, q_1) which results in the :math:`X \otimes Z` tensor order.
        Instead, if we apply it on (q_1, q_0), the matrix will
        be :math:`Z \otimes X`:

        .. parsed-literal::

                 ┌─────────┐
            q_0: ┤1        ├
                 │  Rzy(θ) │
            q_1: ┤0        ├
                 └─────────┘

        .. math::

            \newcommand{\rotationangle}{\frac{\theta}{2}}

            R_{ZY}(\theta)\ q_1, q_0 = exp(-i \frac{\theta}{2} Z{\otimes}Y) =
                \begin{pmatrix}
                    \cos(\rotationangle)   & -\sin(\rotationangle) & 0           & 0          \\
                    \sin(\rotationangle) & \cos(\rotationangle)   & 0           & 0          \\
                    0           & 0           & \cos(\rotationangle)   & \sin(\rotationangle) \\
                    0           & 0           & -\sin(\rotationangle)  & \cos(\rotationangle)
                \end{pmatrix}
    """

    def __init__(
        self,
        theta: ParameterValueType,
        label: Optional[str] = None,
        *,
        duration=None,
        unit="dt"
    ):
        """Create new RZY gate."""
        super().__init__("rzy", 2, [theta], label=label, duration=duration, unit=unit)

    def _define(self):
        """
        gate rzy(theta) a, b { h b; cx a, b; u1(theta) b; cx a, b; h b;}
        """
        # pylint: disable=cyclic-import
        from qiskit.circuit.quantumcircuit import QuantumCircuit
        from qiskit.circuit.library.standard_gates.rx import RXGate
        from qiskit.circuit.library.standard_gates.x import CXGate
        from qiskit.circuit.library.standard_gates.rz import RZGate

        # q_0: ──────────────■─────────────■─────────────
        #      ┌──────────┐┌─┴─┐┌───────┐┌─┴─┐┌─────────┐
        # q_1: ┤ Rx(-π/2) ├┤ X ├┤ Rz(0) ├┤ X ├┤ Rx(π/2) ├
        #      └──────────┘└───┘└───────┘└───┘└─────────┘
        theta = self.params[0]
        q = QuantumRegister(2, "q")
        qc = QuantumCircuit(q, name=self.name)
        rules = [
            (RXGate(-np.pi / 2), [q[1]], []),
            (CXGate(), [q[0], q[1]], []),
            (RZGate(theta), [q[1]], []),
            (CXGate(), [q[0], q[1]], []),
            (RXGate(np.pi / 2), [q[1]], []),
        ]
        for instr, qargs, cargs in rules:
            qc._append(instr, qargs, cargs)

        self.definition = qc

    def inverse(self, annotated: bool = False):
        """Return inverse RZX gate (i.e. with the negative rotation angle).

        Args:
            annotated: when set to ``True``, this is typically used to return an
                :class:`.AnnotatedOperation` with an inverse modifier set instead of a concrete
                :class:`.Gate`. However, for this class this argument is ignored as the inverse
                of this gate is always a :class:`.RZXGate` with an inverted parameter value.

         Returns:
            RZXGate: inverse gate.
        """
        return RZYGate(-self.params[0])

    def __array__(self, dtype=None):
        """Return a numpy.array for the RZX gate."""

        half_theta = float(self.params[0]) / 2
        cos = math.cos(half_theta)
        sin = math.sin(half_theta)
        return np.array(
            [
                [cos, 0, -sin, 0],
                [0, cos, 0, sin],
                [sin, 0, cos, 0],
                [0, -sin, 0, cos],
            ],
            dtype=dtype,
        )

    def power(self, exponent: float):
        """Raise gate to a power."""
        (theta,) = self.params
        return RZYGate(exponent * theta)

    def __eq__(self, other):
        if isinstance(other, RZYGate):
            return self._compare_parameters(other)
        return False

from qiskit import QuantumCircuit
from qiskit.circuit import Parameter

def dig_count_ansatz(pauli_op):
    # Create a quantum circuit with the same number of qubits as the pauli_op
    num_qubits = pauli_op.num_qubits
    qc = QuantumCircuit(num_qubits)

    # Iterate over the terms in the pauli_op
    for label, coeff in zip(pauli_op.paulis.to_labels(), pauli_op.coeffs):
        y_indices = [i for i, x in enumerate(label) if x == 'Y']
        z_indices = [i for i, x in enumerate(label) if x == 'Z']

        # Add parameterized Rz gates for single 'Y' terms
        if len(y_indices) == 1 and len(z_indices) == 0:
            theta = Parameter(f"theta_{y_indices[0]}")
            qc.rz(theta, y_indices[0])

        # Add parameterized RZY gates for 'Y' and 'Z' terms
        elif len(y_indices) == 1 and len(z_indices) == 1:
            theta = Parameter(f"theta_{y_indices[0]}_{z_indices[0]}")
            # Assuming RZY gate accepts parameters and acts on two qubits
            qc.append(RZYGate(theta), [y_indices[0], z_indices[0]])

    return qc


