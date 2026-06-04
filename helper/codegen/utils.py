class Indentable(str):
    """Wrapper class to easily indent code blocks"""
    def __new__(cls, *args, **kwargs):
        return super(cls, cls).__new__(cls, *args, **kwargs)

    def __rshift__(self, level):
        return level * " " + self.replace("\n", "\n" + level * " ")
