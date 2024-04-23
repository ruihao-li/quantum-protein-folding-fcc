"""Gathers penalty parameters for an FCC Hamiltonian."""


class PenaltyParameters:

    def __init__(
        self,
        penalty_back: float = 10.0,
        penalty_redun: float = 10.0,
        penalty_olap_1: float = 10.0,
        penalty_olap_2: float = 5.0,
        penalty_olap_3: float = 10.0,
        penalty_olap_4: float = 5.0,
    ):
        """
        Args:
            penalty_back: A penalty parameter used to penalize consecutive turns that are in the opposite directions.
            penalty_redun: A penalty parameter used to penalize redundant bitstring configurations that do not correspond to any physical turns.
            penalty_olap_1: A penalty parameter used to penalize long-range overlaps using "unbalanced penalization" (arXiv:2211.13914).
            penalty_olap_2: A penalty parameter used to penalize long-range overlaps using "unbalanced penalization" (arXiv:2211.13914).
            penalty_olap_3: A penalty parameter used to penalize long-range overlaps using "unbalanced penalization" (arXiv:2211.13914).
            penalty_olap_4: A penalty parameter used to penalize long-range overlaps using "unbalanced penalization" (arXiv:2211.13914).
        """

        self._penalty_back = penalty_back
        self._penalty_redun = penalty_redun
        self._penalty_olap_1 = penalty_olap_1
        self._penalty_olap_2 = penalty_olap_2
        self._penalty_olap_3 = penalty_olap_3
        self._penalty_olap_4 = penalty_olap_4

    @property
    def penalty_back(self) -> float:
        """Returns a penalty parameter used to penalize consecutive turns that are in the opposite directions."""
        return self._penalty_back

    @property
    def penalty_redun(self) -> float:
        """Returns a penalty parameter used to penalize redundant bitstring configurations that do not correspond to any physical turns."""
        return self._penalty_redun

    @property
    def penalty_olap_1(self) -> float:
        """Returns a penalty parameter used to penalize long-range overlaps using "unbalanced penalization" (arXiv:2211.13914)."""
        return self._penalty_olap_1

    @property
    def penalty_olap_2(self) -> float:
        """Returns a penalty parameter used to penalize long-range overlaps using "unbalanced penalization" (arXiv:2211.13914)."""
        return self._penalty_olap_2

    @property
    def penalty_olap_3(self) -> float:
        """Returns a penalty parameter used to penalize long-range overlaps using "unbalanced penalization" (arXiv:2211.13914)."""
        return self._penalty_olap_3

    @property
    def penalty_olap_4(self) -> float:
        """Returns a penalty parameter used to penalize long-range overlaps using "unbalanced penalization" (arXiv:2211.13914)."""
        return self._penalty_olap_4
