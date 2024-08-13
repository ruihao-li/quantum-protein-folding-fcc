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

import math
import numpy as np
from typing import Optional
from qiskit.circuit.gate import Gate
from qiskit.circuit.quantumregister import QuantumRegister
from qiskit.circuit.parameterexpression import ParameterValueType


class RXYGate(Gate):
    r"""A parametric 2-qubit :math:`X \otimes Y` interaction (rotation about XY).

    **Circuit Symbol:**

    .. parsed-literal::

             ┌─────────┐
        q_0: ┤0        ├
             │  Rxy(θ) │
        q_1: ┤1        ├
             └─────────┘

    **Matrix Representation:**

    .. math::

        R_{XY}(\theta)\ q_0, q_1 = \exp\left(-i \frac{\theta}{2} X{\otimes}Y\right) =
            \begin{pmatrix}
                \cos\left(\theta/2\right) & 0 & 0 & -i\sin\left(\theta/2\right) \\
                0 & \cos\left(\theta/2\right) & i\sin\left(\theta/2\right) & 0 \\
                0 & i\sin\left(\theta/2\right) & \cos\left(\theta/2\right) & 0 \\
                -i\sin\left(\theta/2\right) & 0 & 0 & \cos\left(\theta/2\right)
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
        """Create new RXY gate."""
        super().__init__("rxy", 2, [theta], label=label, duration=duration, unit=unit)

    def _define(self):
        """
        gate rxy(theta) a, b { h a; h b; s b; cx a, b; u1(theta) b; cx a, b; sdg b; h b; h a; }
        """
        from qiskit.circuit.quantumcircuit import QuantumCircuit
        from qiskit.circuit.library.standard_gates.rx import RXGate
        from qiskit.circuit.library.standard_gates.x import CXGate
        from qiskit.circuit.library.standard_gates.s import SGate, SdgGate
        from qiskit.circuit.library.standard_gates.rz import RZGate
        from qiskit.circuit.library.standard_gates.h import HGate

        theta = self.params[0]
        q = QuantumRegister(2, "q")
        qc = QuantumCircuit(q, name=self.name)
        rules = [
            (HGate(), [q[0]], []),
            (HGate(), [q[1]], []),
            (SGate(), [q[1]], []),
            (CXGate(), [q[0], q[1]], []),
            (RZGate(theta), [q[1]], []),
            (CXGate(), [q[0], q[1]], []),
            (SdgGate(), [q[1]], []),
            (HGate(), [q[1]], []),
            (HGate(), [q[0]], []),
        ]
        for instr, qargs, cargs in rules:
            qc._append(instr, qargs, cargs)

        self.definition = qc

    def inverse(self):
        """Return inverse RXY gate (i.e., with the negative rotation angle)."""
        return RXYGate(-self.params[0])

    def __array__(self, dtype=None):
        """Return a numpy.array for the RXY gate."""
        half_theta = float(self.params[0]) / 2
        cos = math.cos(half_theta)
        sin = math.sin(half_theta)
        return np.array([
            [cos, 0, 0, -1j*sin],
            [0, cos, 1j*sin, 0],
            [0, 1j*sin, cos, 0],
            [-1j*sin, 0, 0, cos]
        ], dtype=dtype)

    def power(self, exponent: float):
        """Raise gate to a power."""
        (theta,) = self.params
        return RXYGate(exponent * theta)

    def __eq__(self, other):
        if isinstance(other, RXYGate):
            return self._compare_parameters(other)
        return False



from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.circuit.library import RZZGate, RXXGate, RYYGate, RZXGate
from qiskit.quantum_info import SparsePauliOp
# Import or define other necessary gates, e.g., RZYGate if available

def dig_count_ansatz(pauli_op: SparsePauliOp):
    num_qubits = pauli_op.num_qubits
    qc = QuantumCircuit(num_qubits)
    parameter_dict = {}  # Dictionary to track parameters

    for idx, (label, coeff) in enumerate(zip(pauli_op.paulis.to_labels(), pauli_op.coeffs)):
        x_indices = [i for i, x in enumerate(label) if x == 'X']
        y_indices = [i for i, x in enumerate(label) if x == 'Y']
        z_indices = [i for i, x in enumerate(label) if x == 'Z']

        # Single-qubit gates with unique parameters
        for i in x_indices:
            param_name = f"theta_X_{i}_{idx}"
            if param_name not in parameter_dict:
                parameter_dict[param_name] = Parameter(param_name)
            qc.rx(parameter_dict[param_name], i)
        for i in y_indices:
            param_name = f"theta_Y_{i}_{idx}"
            if param_name not in parameter_dict:
                parameter_dict[param_name] = Parameter(param_name)
            qc.ry(parameter_dict[param_name], i)
        for i in z_indices:
            param_name = f"theta_Z_{i}_{idx}"
            if param_name not in parameter_dict:
                parameter_dict[param_name] = Parameter(param_name)
            qc.rz(parameter_dict[param_name], i)

        # Two-qubit gates with unique parameters
        # Ensuring unique parameter names for each gate application
        if len(x_indices) == 1 and len(y_indices) == 1 and not z_indices:
            param_name = f"theta_RXY_{x_indices[0]}_{y_indices[0]}_{idx}"
            if param_name not in parameter_dict:
                parameter_dict[param_name] = Parameter(param_name)
            qc.append(RXYGate(parameter_dict[param_name]), [x_indices[0], y_indices[0]])
        if len(x_indices) == 1 and len(z_indices) == 1 and not y_indices:
            param_name = f"theta_RXZ_{x_indices[0]}_{z_indices[0]}_{idx}"
            if param_name not in parameter_dict:
                parameter_dict[param_name] = Parameter(param_name)
            qc.append(RZXGate(parameter_dict[param_name]), [x_indices[0], z_indices[0]])
        if len(y_indices) == 1 and len(z_indices) == 1 and not x_indices:
            param_name = f"theta_RYZ_{y_indices[0]}_{z_indices[0]}_{idx}"
            if param_name not in parameter_dict:
                parameter_dict[param_name] = Parameter(param_name)
            qc.append(RZYGate(parameter_dict[param_name]), [y_indices[0], z_indices[0]])

        if len(z_indices) == 2 and not x_indices and not y_indices:
            param_name = f"theta_RZZ_{z_indices[0]}_{z_indices[1]}_{idx}"
            if param_name not in parameter_dict:
                parameter_dict[param_name] = Parameter(param_name)
            qc.append(RZZGate(parameter_dict[param_name]), [z_indices[0], z_indices[1]])
        if len(y_indices) == 2 and not x_indices and not z_indices:
            param_name = f"theta_RYY_{y_indices[0]}_{y_indices[1]}_{idx}"
            if param_name not in parameter_dict:
                parameter_dict[param_name] = Parameter(param_name)
            qc.append(RYYGate(parameter_dict[param_name]), [y_indices[0], y_indices[1]])

    return qc



