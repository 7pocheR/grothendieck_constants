"""Reconstruct the complete rational inventory by two different finite expansions."""
from fractions import Fraction as Q
from package_support import SOURCE, candidate, compressed, load, require, verify_package
from full_common import inventory

if __name__ == '__main__':
    verify_package()
    task = load(SOURCE/'full_task.json')
    d = candidate()
    first = inventory(d, Q(task['threshold']))
    second = inventory(d, Q(task['threshold']), independently=True)
    require(first == second == compressed('inventory.json.gz'), 'Complete inventories differ')
    print('PASS_COMPLETE_RATIONAL_INVENTORY: 120470 modes, 9948 retained, 110522 omissions charged')
