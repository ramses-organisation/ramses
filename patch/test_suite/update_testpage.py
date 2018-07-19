import sys
import glob
from datetime import date, timedelta, datetime


MonthsDict = {'Jan': 1, 'Feb': 2, 'Mar': 3,
              'Apr': 4, 'May': 5, 'Jun': 6,
              'Jul': 7, 'Aug': 8, 'Sep': 9,
              'Oct': 10, 'Nov': 11, 'Dec': 12}

MonthsDictReverse = dict((v, k) for k, v in MonthsDict.iteritems())

def wiki_entry(date, commit, success):
    """Format extracted info as a Markdown table entry"""

    # format wiki entry
    pdf_string = ('[pdf](daily_tests/{year}-{month:02}-{day:02}.pdf)'
                  .format(year=date.year, month=date.month, day=date.day))
    log_string = ('[log](daily_tests/{year}-{month:02}-{day:02}.log)'
                  .format(year=date.year, month=date.month, day=date.day))
    commit_string = ('[{short_hash}](https://bitbucket.org/rteyssie/ramses/commits/{long_hash})'
            .format(short_hash=commit[0:7], long_hash=commit))
    if success:
        success_string = '![ok](ok.png)'
    else:
        success_string = '![fail](fail.png)'

    return str(date)+' '+pdf_string+' '+log_string+' '+commit_string+' '+success_string+' |'


def extract_test_info(filename):
    """Extracts crucial test info from the log file"""

    with open(filename) as infile:
        passed = True
        last_commit = 7*'?'

        # figure out the date from the filename
        compile_date = date(int(filename[-14:-10]),
                            int(filename[-9:-7]),
                            int(filename[-6:-4]))

        for line in infile:
            if 'compile' in line:
                    passed = True
            # find if commit info exists
            if 'commit' in line.split(' '):
                last_commit = line.split('=')[-1][1:-1]

            # if any of the tests failed, then all failed
            if 'failed' in line.split(' '):
                passed = False

    return compile_date, last_commit, passed


def rebuild_wiki(logs, wiki, max_months=12):
    """Rebuilts the AutoTests.md wiki file from scratch"""

    row = ''
    prev_month = 99
    num_months = 0
    dates = []
    commits = []
    successes = []


    months = []
    commits2 = []
    successes2 = []

    # process logs
    for logname in logs:
        compile_date, last_commit, passed = extract_test_info(logname)

        dates.append(compile_date)
        commits.append(last_commit)
        successes.append(passed)

    today = datetime.today()
    today = date(today.year, today.month, today.day)
    x=len(dates)
    for i in xrange(x):
        current = dates[~i]
        if (today - current).days < 365:
            if prev_month != current.month:
                if i == 0:
                    months.append(dates[~i-current.day+1:])
                    commits2.append(commits[~i-current.day+1:])
                    successes2.append(successes[~i-current.day+1:])
                    prev_month = current.month

                else:
                    months.append(dates[~i-current.day+1:~i+1])
                    commits2.append(commits[~i-current.day+1:~i+1])
                    successes2.append(successes[~i-current.day+1:~i+1])
                    prev_month = current.month

    # buld a new wiki md file
    with open(wiki, 'w') as wikifile:
        wikifile.write('# RAMSES Daily Test\n')
        wikifile.write('\n')
        wikifile.write('This page contains the test results for the last year from the RAMSES test suite, which is run daily. The test pipeline was provided by Neil Vaytet with fixes/expansions by Pawel Biernacki.\n')
        wikifile.write('\n')
        # table entries

        for j in xrange(len(months)):
            new_month = True
            for i in xrange(len(months[j])):
                current = months[j][i]
                last = months[j][-1]

                if months[j][i].month != last.month:
                    continue

                if new_month:
                    wikifile.write('# %s %d' % (MonthsDictReverse[last.month], last.year))
                    wikifile.write('\n\n')
                    wikifile.write('| Mon | Tue | Wed | Thu | Fri | Sat | Sun |\n')
                    wikifile.write('| ---:| ---:| ---:| ---:| ---:| ---:| ---:|\n')
                    row = ''
                    for _ in xrange(current.weekday()+1):
                        row += '| '

                row += wiki_entry(months[j][i], commits2[j][i], successes2[j][i])

                # taking care of missing data
                if i+1 < len(months[j]):
                    day_diff = months[j][i+1].day-current.day
                    if day_diff > 1:
                        if day_diff-1 >= 6-current.weekday():
                            row += (6-current.weekday()-1)*'| '
                            row += '\n'
                            row += (day_diff-(6-current.weekday()))*'| '
                        else:
                            row += (day_diff-1)*'| '

                if current.weekday() == 6 or i+1 == len(months[j]):
                    row += '\n'
                    wikifile.write(row)
                    row = ''

                if i+1 == len(months[j]):
                    wikifile.write(3*'\n')

                new_month = False

    return



def add_last_entry(logfile, wiki):
    """Finds last logfile, extracts the test info and adds it to the wiki"""

    # gather info on the last test
    when, commit, passed = extract_test_info(logfile)

    # read the content of wiki
    with open(wiki, 'r') as wikifile:
        wiki_contents = wikifile.readlines()

    # insert the latest test after the header
    wiki_contents.insert(6, wiki_entry(when, commit, passed))

    # write the updated content of the wiki
    with open(wiki, 'w') as wikifile:
        wiki_contents = ''.join(wiki_contents)
        wikifile.write(wiki_contents)


def main():

    wikidir = sys.argv[1]
    wikifile = sys.argv[2]

    # gather all the log files
    logs =  glob.glob(wikidir+'daily_tests/201*.log')

    rebuild_wiki(logs, wikidir+'/'+wikifile)


if __name__ == '__main__':
    main()
