class Run:
    def __init__(self, data):
        self.iteration = data[0]
        self.s = data[1]
        self.x = data[2]
        self.their_value = data[3]
        self.my_value = data[4]
        self.diff = data[5]
        self.error = data[6]

tests = []
data = []
f = open("Log.2.6.2017-14:2:21")

for line in f:
    words = line.split()
    if(len(words) > 0):
        if words[0] == "Iteration":
            data.append(int(words[1][0]))
        if words[0] == "Inputs:":
            data.append(float(words[3]))
            data.append(float(words[6]))
        if words[0] == "Their":
            data.append(float(words[2]))
        if words[0] == "My":
            data.append(float(words[2]))
        if words[0] == "Difference:":
            data.append(float(words[1]))
        if words[0] == "Error:":
            data.append(float(words[1]))
            tests.append(Run(data))
            data = []

errors = []

for run in tests:
    errors.append(run.error)
errors.sort()
