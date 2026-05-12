"""Assignment helper modules for XpongeCPP."""

from collections import OrderedDict


class AssignRule:
    """Xponge-compatible atom typing rule registry.

    Built-in typing rules stay in the C++ core. This class is the compatibility
    path for user-defined Python rules, matching Xponge's add_rule/pre_action/
    post_action/pure_string behavior.
    """

    all = {}

    def __init__(self, name, pure_string=False, pre_action=None, post_action=None):
        self.name = name
        self.rules = OrderedDict()
        self.priority = {}
        self.built = False
        self.pure_string = pure_string
        self.pre_action = pre_action
        self.post_action = post_action
        AssignRule.all[name] = self

    def add_rule(self, atomtype, priority=0):
        def wrapper(rule_function):
            self.rules[str(atomtype)] = rule_function
            self.priority[str(atomtype)] = -priority
            self.built = False
            return rule_function

        return wrapper

    def set_pre_action(self, function):
        self.pre_action = function

    def set_post_action(self, function):
        self.post_action = function


__all__ = ["AssignRule"]
