"""Version-related DGAME experiment constants."""

from dataclasses import dataclass

from dgame.constants import DECKE_LABEL, DIRECTOR_LABEL, MATCHER_LABEL


@dataclass
class DGameVersion:
    version: str
    director_label: str
    non_director_label: str
    n_participant_streams: int

    def participant_roles(self) -> tuple[str, str]:
        return self.director_label, self.non_director_label


SUPPORTED_DGAME_VERSIONS = {
    "2": DGameVersion(
        version="2",
        director_label=DIRECTOR_LABEL,
        non_director_label=DECKE_LABEL,
        n_participant_streams=1,
    ),
    "3": DGameVersion(
        version="3",
        director_label=DIRECTOR_LABEL,
        non_director_label=MATCHER_LABEL,
        n_participant_streams=2,
    ),
}
