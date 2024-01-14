from dataclasses import dataclass


@dataclass
class Student:
    subject: str = "Maths"
    name: str = "Test"


st = Student("Biology", "Jack")

print(st.name, st.subject)
