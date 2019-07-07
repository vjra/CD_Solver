import webbrowser
import os
import time

switch = True
print(len(os.listdir('.')))
while switch == True:

    if len(os.listdir('.')) >2:
        webbrowser.open('https://www.youtube.com/watch?v=BIKarAqOB9I')
        switch = False
        print(switch)
    else:
        print(switch)

    time.sleep(60)
