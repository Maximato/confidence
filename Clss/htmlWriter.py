class htmlWriter:
    def __init__(self, html):
        self.html = html

    def write_to(self, filename):
        with open(filename, "w") as f:
            f.write(self.html)
