with open("debug.txt", "a") as file:
    births = 1
    deaths = 1
    line = str(births) + "," + str(deaths) + "\n"
    file.write(line)
    file.write("hello", 'a')