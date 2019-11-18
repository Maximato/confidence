class HTML:

    classes = {
        "c00": '<span class="c00">',
        "c01": '<span class="c01">',
        "c10": '<span class="c10">',
        "c20": '<span class="c20">',
        "c30": '<span class="c30">',
        "c40": '<span class="c40">',
        "c50": '<span class="c50">',
        "c60": '<span class="c60">',
        "c70": '<span class="c70">',
        "c80": '<span class="c80">',
        "c90": '<span class="c90">'
    }

    @staticmethod
    def get_header(file_name="./header.txt"):
        with open(file_name, "r") as f:
            header = f.read()
        return header

    def create_consensus(self, consensus):
        html = self.get_header()

        symbols = consensus["symbols"]
        groups = consensus["groups"]
        for i in range(len(symbols)):
            html += self.classes[groups[i]] + symbols[i] + "</span>"

        html += "</body>\n</html>"

        return html
    
print(HTML().get_header())
