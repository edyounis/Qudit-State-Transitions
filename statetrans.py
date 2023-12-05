"""
Qudit State Transition Library Code

Designed to find scalable solutions to the qudit GHZ state synthesis problem
without the use of parameterized gates.
"""


import heapq
from typing import Tuple
from typing import Callable


State = Tuple[str, ...]
"""
The type of a state in our system is an unparameterized quantum state.

This is given as a sum of computational basis states, represented as strings.

For example, the three-qutrit ghz state |000> + |111> + |222> would be
represented as ('000', '111', '222').
"""


Transition = Callable[[State], State]
"""
A transition from one state to another marks a valid operation.

Transitions build our state space. We apply transitions during search to
find a solution state. These will typically represent a gate or a sequence
of gates.

For example, the application of an X gate in the subspace spanned by the
|0> and |1> states to the first qudit would transition the state |0XY> to
|1XY>.
"""

Path = Tuple[str, ...]
"""A path is a sequence of transition names."""


class heapnode:
    """Node in Dijkstra's search graph."""

    def __init__(self, state, path, cost):
        self.state = state
        self.path = path
        self.cost = cost

    def __lt__(self, other):
        return self.cost < other.cost

    def __eq__(self, other):
        return self.cost == other.cost

    def __gt__(self, other):
        return self.cost > other.cost


def state_search_via_dijkstras(
    init_state: State,
    transitions: dict[str, Transition],
    success: Callable[[State], bool],
    cost: Callable[[Path], int],
) -> Path:
    """
    Performs search from init to success using Dijkstra's algorithm.

    Args:
        init_state (State): The initial state.

        transitions (dict[str, Transition]): A dictionary from transition
            names to transitions.

        success (Callable[[State], bool]): A function that returns true
            if a state is a success state.

        cost (Callable[[Path], int]): A function that returns the
            cost of a path.
    
    Returns:
        Path: A path from init to success.
    """
    visited = {}
    frontier = []
    heapq.heappush(frontier, heapnode(init_state, tuple(), cost(tuple())))

    while len(frontier) > 0:
        node = heapq.heappop(frontier)

        if node.state in visited:
            continue

        visited[node.state] = node.path

        if success(node.state):
            return node.path
        
        for name, transition in transitions.items():
            new_state = transition(node.state)
            if new_state not in visited:
                new_path = node.path + (name,)
                heapq.heappush(frontier, heapnode(new_state, new_path, cost(new_path)))
    
    raise RuntimeError("No solution found.")


__all__ = ['state_search_via_dijkstras', 'State', 'Transition', 'Path']
