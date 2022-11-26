# This function will take in start and end dates and yield an iteration
# of dates from start to end

# 15July2021


import datetime


def func_daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days) + 1):
        yield start_date + datetime.timedelta(n)
