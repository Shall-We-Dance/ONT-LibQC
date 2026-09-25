"""Public, bounded settings for a physical P5-to-P7 library layout."""
from typing import Literal
from pydantic import Field, model_validator
from pydantic import BaseModel, ConfigDict

P5 = 'AATGATACGGCGACCACCGAGATCTACAC'
P7_FORWARD = 'ATCTCGTATGCCGTCTTCTGCTTG'
POSITIONS = ['before_i5', 'after_i5', 'before_insert', 'after_insert', 'before_i7', 'after_i7']

class UMI(BaseModel):
    model_config = ConfigDict(extra='forbid')
    position: Literal['before_i5', 'after_i5', 'before_insert', 'after_insert', 'before_i7', 'after_i7']
    length: int = Field(strict=True, ge=1, le=64)

class LibraryStructure(BaseModel):
    model_config = ConfigDict(extra='forbid')
    mode: Literal['auto', 'custom'] = 'auto'
    adapter_family: Literal['auto', 'truseq', 'nextera'] = 'auto'
    p5_sequence: str | None = Field(default=None, pattern=r'^[ACGT]{16,80}$')
    p7_sequence: str | None = Field(default=None, pattern=r'^[ACGT]{16,80}$')
    adapter_i5_sequence: str | None = Field(default=None, pattern=r'^[ACGT]{16,80}$')
    adapter_i7_sequence: str | None = Field(default=None, pattern=r'^[ACGT]{16,80}$')
    i5_length: int | None = Field(default=None, strict=True, ge=0, le=24)
    i7_length: int | None = Field(default=None, strict=True, ge=0, le=24)
    umis: list[UMI] = Field(default_factory=list, max_length=6)

    @model_validator(mode='after')
    def valid_layout(self):
        positions = [u.position for u in self.umis]
        if len(set(positions)) != len(positions):
            raise ValueError('Use at most one UMI at each position')
        if self.mode == 'auto' and (self.umis or self.adapter_family != 'auto' or any(v is not None for v in (self.p5_sequence, self.p7_sequence, self.adapter_i5_sequence, self.adapter_i7_sequence, self.i5_length, self.i7_length))):
            raise ValueError('Choose Customize to set a library structure or add UMIs')
        self.umis.sort(key=lambda umi: POSITIONS.index(umi.position))
        return self


def structure_catalog():
    return {'p5_sequence': P5, 'p7_sequence': P7_FORWARD, 'umi_positions': POSITIONS,
            'umi_max_length': 64, 'index_max_length': 24,
            'adapter_sequence_min_length': 16, 'adapter_sequence_max_length': 80}
