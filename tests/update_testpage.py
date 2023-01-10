import sys
import glob
from datetime import date, timedelta, datetime


MonthsDict = {'Jan': 1, 'Feb': 2, 'Mar': 3,
              'Apr': 4, 'May': 5, 'Jun': 6,
              'Jul': 7, 'Aug': 8, 'Sep': 9,
              'Oct': 10, 'Nov': 11, 'Dec': 12}

MonthsDictReverse = dict((v, k) for k, v in MonthsDict.items())


def wiki_entry(date, commit, success):
    """Format extracted info as a Markdown table entry"""

    commit_url = 'https://bitbucket.org/rteyssie/ramses/commits/'
    # format wiki entry
    pdf_str = ('[pdf](daily_tests/{year}-{month:02}-{day:02}.pdf)'
               .format(year=date.year, month=date.month, day=date.day))
    log_str = ('[log](daily_tests/{year}-{month:02}-{day:02}.log)'
               .format(year=date.year, month=date.month, day=date.day))
    commit_str = ('[{short_hash}]({commit_url}{long_hash})'
                  .format(commit_url=commit_url,
                          short_hash=commit[0:7],
                          long_hash=commit))
    if success:
        success_str = '![ok](ok.png)'
    else:
        success_str = '![fail](fail.png)'

    return str(date)+' '+pdf_str+' '+log_str+' '+commit_str+' '+success_str+' |'


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


def rebuild_wiki(logs, wiki):
    """Rebuilts the AutoTests.md wiki file from scratch"""

    row = ''
    prev_month = 99
    all_dates = []
    all_commits = []
    all_successes = []

    months = []
    month_commits = []
    month_successes = []

    # process logs
    for logname in logs:
        compile_date, last_commit, passed = extract_test_info(logname)

        all_dates.append(compile_date)
        all_commits.append(last_commit)
        all_successes.append(passed)

    today = datetime.today()
    today = date(today.year, today.month, today.day)

    for i in range(len(all_dates)):
        # getting the current day of the wiki
        current = all_dates[~i]
        # setting the limit
        if (today - current).days < 365:
            if prev_month != current.month:
                if i == 0:  # edge case
                    months.append(all_dates[~i-current.day+1:])
                    month_commits.append(all_commits[~i-current.day+1:])
                    month_successes.append(all_successes[~i-current.day+1:])
                else:
                    months.append(all_dates[~i-current.day+1:~i+1])
                    month_commits.append(all_commits[~i-current.day+1:~i+1])
                    month_successes.append(all_successes[~i-current.day+1:~i+1])

                prev_month = current.month

    del all_dates, all_commits, all_successes

    # buld a new wiki md file
    with open(wiki, 'w') as wikifile:
        wikifile.write('# RAMSES Daily Test\n')
        wikifile.write('\n')
        wikifile.write('This page contains the test results for the last year from the RAMSES test suite, which is run daily. The test pipeline was provided by Neil Vaytet with fixes/expansions by Pawel Biernacki.\n')
        wikifile.write('\n')
        # table entries

        for j in range(len(months)):
            new_month = True
            for i in range(len(months[j])):
                current = months[j][i]
                last = months[j][-1]

                if months[j][i].month != last.month:
                    continue

                # set the header of the table for the month
                if new_month:
                    wikifile.write('# {} {:d}'
                                   .format(MonthsDictReverse[last.month],
                                           last.year))
                    wikifile.write('\n\n')
                    wikifile.write('| Mon | Tue | Wed | Thu | Fri | Sat | Sun |\n')
                    wikifile.write('|'+7*' ---: |'+'\n')
                    row = ''
                    for _ in range(current.weekday()+1):
                        row += '| '

                row += wiki_entry(months[j][i],
                                  month_commits[j][i],
                                  month_successes[j][i])

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

                # if Sunday or end of the month
                if current.weekday() == 6 or i+1 == len(months[j]):
                    row += '\n'
                    wikifile.write(row)
                    row = ''

                # add spacing between months
                if i+1 == len(months[j]):
                    wikifile.write(3*'\n')

                new_month = False

    return


def main():

    wikidir = sys.argv[1]
    wikifile = sys.argv[2]

    # gather all the log files
    logs = sorted(glob.glob(wikidir+'/daily_tests/20*.log'))

    rebuild_wiki(logs, wikidir+'/'+wikifile)


if __name__ == '__main__':
    main()
