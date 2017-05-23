'''
SSEA: Sample Set Enrichment Analysis
'''

__author__ = "Matthew Iyer, Yashar Niknafs"
__copyright__ = "Copyright 2012-2017"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


class DTColumn:
    def __init__(self, name, searchable=False, orderable=False):
        self.name = name
        self.orderable = orderable
        self.searchable = searchable
        self.search_value = ''


class DTQuery:
    def __init__(self, params, columns):
        self.columns = columns
        self.draw = int(params.get('draw'))
        self.start = int(params.get('start'))
        self.length = int(params.get('length'))
        self.global_search = params.get('search[value]', '')
        self.sort_expressions = []

        # get column params
        for i in range(len(columns)):
            filter_expr = None
            searchable = params.get('columns[{:d}][searchable]'.format(i), '')
            orderable = params.get('columns[{:d}][orderable]'.format(i), '')
            search_value = params.get('columns[{:d}][search][value]'.format(i), '')
            self.columns[i].searchable = (searchable == 'true')
            self.columns[i].orderable = (orderable == 'true')
            self.columns[i].search_value = search_value

        # order by columns
        i = 0
        while params.get('order[{:d}][column]'.format(i), False):
            column_index = int(params.get('order[{:d}][column]'.format(i)))
            direction = params.get('order[{:d}][dir]'.format(i))
            column = columns[column_index]
            if column.name:
                self.sort_expressions.append((column.name, direction))
            i += 1

    def execute(self, conn, table):
        c = conn.cursor()

        # count total records
        c.execute('select count(*) from {t}'.format(t=table))
        records_total = c.fetchone()[0]

        filter_expressions = []
        # global search
        if self.global_search:
            val = "'%{}%'".format(self.global_search)
            filter_expr = []
            for column in self.columns:
                if column.searchable:
                    filter_expr.append('{} {} {}'.format(column.name, 'like', val))
            filter_expr = '({})'.format(' or '.join(filter_expr))
            filter_expressions.append(filter_expr)

        # column search
        for column in self.columns:
            if column.search_value:
                val = "%{}%".format(column.search_value)
                filter_expr = '{} {} {}'.format(column.name, 'like', val)
                filter_expressions.append(filter_expr)

        # filter string
        where_str = ''
        records_filtered = records_total
        if len(filter_expressions) > 0:
            where_str = ' where {}'.format(' and '.join(filter_expressions))
            # compute number of filtered records
            c.execute('select count(*) from {t} {where}'.format(t=table, where=where_str))
            records_filtered = c.fetchone()[0]

        # sort expression
        orderby_str = ''
        if len(self.sort_expressions) > 0:
            orderby_str = ' order by {}'.format(','.join('{} {}'.format(*expr) for expr in self.sort_expressions))

        # query
        stmt = 'select {columns} from {table}{where}{orderby} limit {limit} offset {offset}'
        columns_str = ','.join(column.name for column in self.columns if column.name)
        stmt = stmt.format(columns=columns_str,
                           table=table,
                           where=where_str,
                           orderby=orderby_str,
                           limit=self.length,
                           offset=self.start)
        c = conn.execute(stmt)
        entries = c.fetchall()
        c.close()

        # setup output
        error = False
        output = {}
        output['draw'] = str(self.draw)
        output['recordsTotal'] = records_total
        output['recordsFiltered'] = records_filtered
        if error:
            output['error'] = error
            return output
        column_names = [column.name for column in self.columns if column.name]
        output['data'] = [{k: v for k, v in zip(column_names, row)}
                          for row in entries]
        return output
