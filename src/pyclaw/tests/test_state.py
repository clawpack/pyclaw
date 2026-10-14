import pytest
from pyclaw import State
import pyclaw
from clawpack.pyclaw import Patch
import numpy as np
from examples.advection_2d_annulus import advection_annulus

class TestState():

    @pytest.fixture
    def state(self) -> State:
        r_lower = 0.2
        r_upper = 1.0
        m_r = 40

        theta_lower = 0.0
        theta_upper = np.pi*2.0
        m_theta = 120

        r     = pyclaw.Dimension(r_lower,r_upper,m_r,name='r')
        theta = pyclaw.Dimension(theta_lower,theta_upper,m_theta,name='theta')
        
        patch = Patch([r, theta])
        patch.grid.mapc2p = advection_annulus.mapc2p_annulus
        patch.grid.num_ghost = 2
        num_eqn = 1
        state = State(patch,num_eqn)

        advection_annulus.qinit(state)

        dx, dy = state.grid.delta
        p_corners = state.grid.p_nodes
        state.aux = advection_annulus.edge_velocities_and_area(p_corners[0],p_corners[1],dx,dy)
        state.index_capa = 2
        return state


    @pytest.mark.parametrize(("raise_error", "aux_none", "invalid_index_capa"), 
                             [(False, False, -1), # Valid state
                              (True, True, -1), # Invalid state: aux is None
                              (True, False, -1), # Invalid state: index_capa is -1
                              (True, False, -10)]) # Invalid state: index_capa is out of bounds
    def test_state_initialization(self, state: State, 
                                  raise_error: bool, 
                                  aux_none: bool, 
                                  invalid_index_capa: int) -> None:
        if raise_error:
            if aux_none:
                state.aux = None
            elif state.aux is not None:
                state.index_capa = invalid_index_capa
            with pytest.raises(ValueError):
                assert not state.is_valid()

        else:
            assert state.is_valid()