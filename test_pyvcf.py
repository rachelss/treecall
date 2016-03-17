import pytest
from treecall.pyvcf import *

def test_PL10_order():
    o = PL10_order('A', 'C')
    assert o == [0, 1, 3, 3, 2, 4, 4, 5, 5, 5]
    
    o = PL10_order('A', 'C,T')
    assert o == [0, 1, 6, 3, 2, 7, 4, 9, 6, 5]
    
def test_DPR4_order():
    o = DPR4_order('A', 'C')
    assert o == [0, 1, 2, 2]
    
    o = DPR4_order('A', 'C,T')
    assert o == [0, 1, 3, 2]