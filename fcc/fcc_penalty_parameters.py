"""Gathers penalty parameters for an FCC Hamiltonian."""


class PenaltyParameters:

    def __init__(
        self,
        penalty_back: float = 10.0,
        penalty_redun: float = 10.0,
        penalty_olap: float = 10.0,
    ):
        """
        Args:
            penalty_back: A penalty parameter used to penalize consecutive turns that are in the opposite directions.
            penalty_redun: A penalty parameter used to penalize redundant bitstring configurations that do not correspond to any physical turns.
            penalty_olap: A penalty parameter used to penalize long-range overlaps using higher-order polynomials.
        """

        self._penalty_back = penalty_back
        self._penalty_redun = penalty_redun
        self._penalty_olap = penalty_olap

    @property
    def penalty_back(self) -> float:
        """Returns a penalty parameter used to penalize consecutive turns that are in the opposite directions."""
        return self._penalty_back

    @property
    def penalty_redun(self) -> float:
        """Returns a penalty parameter used to penalize redundant bitstring configurations that do not correspond to any physical turns."""
        return self._penalty_redun

    @property
    def penalty_olap(self) -> float:
        """Returns a penalty parameter used to penalize long-range overlaps using higher-order polynomials."""
        return self._penalty_olap
