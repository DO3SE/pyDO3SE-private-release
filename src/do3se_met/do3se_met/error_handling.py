class InputError(ValueError):

    def __init__(self, field, message):
        self.field = field
        self.message = message

    def __str__(self):
        return f"""Error occured with input field: {self.field}

        !!--------------!!
        {str(self.message)}

        """
