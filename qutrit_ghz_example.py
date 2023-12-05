"""Script to solve for a three-qutrit kernel of the Qutrit GHZ state."""

# ###################################
# Transitions for Qutrit GHZ state
# ###################################

from functools import partial
from statetrans import State

def x(state: State, qubit: int, low_radix: int):
    """
    Apply a subspace X gate to `qubit` in the `state`.

    The subspace is determined by the `low_radix` parameter. The subspace is
    `low_radix` and `low_radix + 1`. For example, if `low_radix` is 0, then the
    subspace is 0 and 1.
    """
    new_state = []

    high_radix = low_radix + 1
    low_radix_str = str(low_radix)
    high_radix_str = str(high_radix)

    for string in state:
        new_string = string[:qubit]

        if string[qubit] == high_radix_str:
            new_string += low_radix_str

        elif string[qubit] == low_radix_str:
            new_string += high_radix_str

        else:
            new_string += string[qubit]

        new_string += string[qubit+1:]
        new_state.append(new_string)
    return tuple(new_state)

def subswap(state: State, low_qubit: int, low_radix: int):
    """
    Apply a subspace SWAP gate to `low_qubit` and `low_qubit + 1` in the `state`.

    The subspace is determined by the `low_radix` parameter. The subspace is
    `low_radix` and `low_radix + 1`. For example, if `low_radix` is 0, then the
    subspace is 0 and 1.
    """
    new_state = []
    high_radix = low_radix + 1
    low_radix_str = str(low_radix) + str(low_radix)
    high_radix_str = str(high_radix) + str(high_radix)

    for string in state:
        new_string = string[:low_qubit]
        if string[low_qubit:low_qubit+2] == low_radix_str:
            new_string += high_radix_str
        elif string[low_qubit:low_qubit+2] == high_radix_str:
            new_string += low_radix_str
        else:
            new_string += string[low_qubit:low_qubit+2]

        new_string += string[low_qubit+2:]
        new_state.append(new_string)
    return tuple(new_state)


x0_01 = partial(x, qubit=0, low_radix=0)
x1_01 = partial(x, qubit=1, low_radix=0)
x2_01 = partial(x, qubit=2, low_radix=0)
x0_12 = partial(x, qubit=0, low_radix=1)
x1_12 = partial(x, qubit=1, low_radix=1)
x2_12 = partial(x, qubit=2, low_radix=1)
subswap12_01 = partial(subswap, low_qubit=1, low_radix=0)
subswap12_12 = partial(subswap, low_qubit=1, low_radix=1)


transitions = {
    "x1_01": x1_01,
    "x2_01": x2_01,
    "x1_12": x1_12,
    "x2_12": x2_12,
    "subswap12_01": subswap12_01,
    "subswap12_12": subswap12_12,
}


# #############################
# Search for Qutrit GHZ state
# #############################

from statetrans import state_search_via_dijkstras


success = lambda state: '000' in state and '111' in state and '222' in state
"""
The success state is the Qutrit GHZ state, represented as ('000', '111', '222').
We don't care about the order of the states, so we just check for inclusion.
Our transitions will also not change the number of states, so we don't need to
check for exclusion of other states.
"""


init_state = ('000', '110', '220')
"""
The initial state is the three-qutrit bell state |000> + |110> + |220>.

Once we solve for the Qutrit GHZ state, we can use the same transitions to
build a Qutrit GHZ state of any size. This is because we can prepare the
bell state in first two qutrits, and then use the transitions to create a
three-qutrit GHZ state. We can cycle out the first qutrit and in the next
qutrit, which will also be in be in the bell state, and then use the
transitions to create a four-qutrit GHZ state. We can repeat this process
until we have a GHZ state of any size. This, of course, ignores leading
phases.
"""


cost = lambda path: len([x for x in path if x.startswith('subswap')]) * 10 + len(path)
"""
The cost is dominated by the number of two-qudit gates, but we also add a
penalty for the number of one-qudit gates.
"""


solution = state_search_via_dijkstras(init_state, transitions, success, cost)


for transition in solution:
    print(transition)

