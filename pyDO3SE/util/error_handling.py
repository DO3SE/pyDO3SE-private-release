
class ConfigError(Exception):
    pass


class DayRangeError(ConfigError):
    pass


class InputDataError(Exception):
    pass


class OutputError(Exception):
    pass

class Do3seRunError(Exception):
    def __init__(self, message, error=None):
        super().__init__(message)
        self.error = error


class Do3seSetupError(Exception):
    pass