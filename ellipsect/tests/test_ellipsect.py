

import pytest

from ellipsect import ellipsectors 


# simple run test
def test_exit():
    with pytest.raises(SystemExit) as e:
        ellipsectors.run()
    assert e.type == SystemExit 
    assert e.value.code == 0 


