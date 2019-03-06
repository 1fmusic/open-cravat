#!/usr/bin/env python3

import argparse
import os
import sys
import yaml
import sqlite3
import json
import re
import time

def get_keycol (level):
    if level == 'variant':
        keycol = 'base__uid'
    elif level == 'gene':
        keycol = 'base__hugo'
    return keycol

class FilterColumn(object):

    test2sql = {
        'equals': '==',
        'lessThanEq': '<=',
        'lessThan': '<',
        'greaterThanEq': '>=',
        'greaterThan': '>',
        'hasData': 'is not null',
        'noData': 'is null',
        'stringContains': 'like',
        'stringStarts': 'like',
        'stringEnds': 'like',
        'between': 'between',
        'in': 'in',
        'select': 'in',
    }

    def __init__(self, d, parent_operator):
        self.column = d['column']
        self.test = d['test']
        self.value = d.get('value')
        self.negate = d.get('negate', False)
        self.parent_operator = parent_operator

    def get_sql(self):
        incexc = 'include'
        if self.column == 'tagsampler__samples' and type(self.value) == list:
            if type(self.value) == list:
                s = 's.base__sample_id="' + self.value[0] + '"'
                for v in self.value[1:]:
                    s += ' or s.base__sample_id="' + v + '"'
            elif type(self.value) == str:
                s = 's.base__sample_id="' + self.value + '"'
            if self.negate and self.parent_operator == 'AND':
                incexc = 'exclude'
        elif self.column == 'tagsampler__tags':
            s = 'm.base__tags="' + self.value[0] + '"'
            for v in self.value[1:]:
                s += ' or m.base__tags="' + v + '"'
        elif self.test == 'multicategory':
            s = 't.{} like "%{}%"'.format(self.column, self.value[0])
            for v in self.value[1:]:
                s += ' or t.{} like "%{}%"'.format(self.column, v)
        else:
            s = 't.{col} {opr}'.format(col=self.column, opr=self.test2sql[self.test])
            sql_val = None
            if self.test == 'equals':
                if type(self.value) is list:
                    v = self.value[0]
                    if type(v) is str:
                        sql_val = '"' + v + '"'
                    else:
                        sql_val = str(v)
                    for v in self.value[1:]:
                        if type(v) is str:
                            v = '"' + v + '"'
                        else:
                            v = str(v)
                        sql_val += ' OR {} == {}'.format(self.column, v)
                else:
                    if type(self.value) is str:
                        sql_val = '"{}"'.format(self.value)
                    else:
                        sql_val = str(self.value)
            elif self.test == 'stringContains':
                sql_val = '"%{}%"'.format(self.value)
            elif self.test == 'stringStarts':
                sql_val = '"{}%"'.format(self.value)
            elif self.test == 'stringEnds':
                sql_val = '"%{}"'.format(self.value)
            elif self.test in ('select','in'):
                str_toks = []
                for val in self.value:
                    if type(val) == str:
                        str_toks.append('"{}"'.format(val))
                    else:
                        str_toks.append(str(val))
                sql_val = '(' + ', '.join(str_toks) + ')'
            elif self.test == 'between':
                sql_val = '{} and {}'.format(self.value[0], self.value[1])
            elif self.test in ('lessThan','lessThanEq','greaterThan','greaterThanEq'):
                sql_val = str(self.value)
            if sql_val:
                s += ' '+sql_val
        if self.negate and incexc != 'exclude':
            s = 'not('+s+')'
        return s, incexc

class FilterGroup(object):
    def __init__(self, d, level):
        self.level = level
        self.operator = d.get('operator', 'and')
        self.negate = d.get('negate',False)
        self.groups = [FilterGroup(x, self.level) for x in d.get('groups',[])]
        self.columns = [FilterColumn(x, self.operator) for x in d.get('columns', [])]

    def determine_sample_or_tag_needed (self, sqls):
        sample_needed = False
        for sql in sqls:
            if 's.base__sample_id' in sql:
                sample_needed = True
                break
        tag_needed = False
        for sql in sqls:
            if 'm.base__tags' in sql:
                tag_needed = True
                break
        from_add = ''
        if sample_needed:
            from_add += ', sample as s'
        if tag_needed:
            from_add += ', mapping as m'
        where_add = ''
        if sample_needed:
            where_add += ' and s.base__uid=t.base__uid'
        if tag_needed:
            where_add += ' and m.base__uid=t.base__uid'
        return from_add, where_add

    def get_sql(self):
        column_include_sqls = []
        column_exclude_sqls = []
        for operand in self.columns:
            sql, incexc = operand.get_sql()
            if sql == '':
                continue
            if incexc == 'include':
                column_include_sqls.append(sql)
            elif incexc == 'exclude':
                column_exclude_sqls.append(sql)
        print(column_include_sqls)
        print(column_exclude_sqls)
        q = ''
        keycol = get_keycol(self.level)
        if len(column_include_sqls) > 0 or len(column_exclude_sqls) > 0:
            q = 'select distinct(t.{}) from {} as t'.format(keycol, self.level)
            if len(column_include_sqls) > 0:
                from_add, where_add = self.determine_sample_or_tag_needed(column_include_sqls)
                q += from_add
                s = ''
                sql_operator = ' ' + self.operator + ' '
                s += sql_operator.join([sql for sql in column_include_sqls])
                q += ' where ({})'.format(s)
                q += where_add
            if len(column_exclude_sqls) > 0:
                from_add, where_add = self.determine_sample_or_tag_needed(column_exclude_sqls)
                q += ' except select distinct(t.{}) from {} as t'.format(keycol, self.level)
                q += from_add
                s = ''
                sql_operator = ' ' + self.operator + ' '
                s += sql_operator.join([sql for sql in column_exclude_sqls])
                q += ' where ({})'.format(s)
                q += where_add
        print('after columns q=', q)
        group_qs = []
        for operand in self.groups:
            group_q = operand.get_sql()
            group_qs.append(group_q)
        if len(group_qs) > 0:
            for group_q in group_qs:
                if q != '':
                    if self.operator == 'AND':
                        q += ' INTERSECT'
                    elif self.operator == 'OR':
                        q += ' UNION'
                q += ' select * from ({})'.format(group_q)
        print('negate=', self.negate)
        if self.negate:
            print('negating. q=', q)
            q = 'select t.{} from {} as t except {}'.format(keycol, self.level, q)
            print('after negating q=', q)
        print('final q=', q)
        return q

class CravatFilter ():
    def __init__ (self, dbpath=None, filterpath=None, filtername=None, 
            filterstring=None, filter=None, mode='sub'):
        self.mode = mode
        if self.mode == 'main':
            self.stdout = True
        else:
            self.stdout = False
        self.dbpath = None
        self.filterpath = None
        self.cmd = None
        self.level = None
        self.filter = None
        self.filtertable = 'filter'
        self.savefiltername = None
        self.filtername = None
        self.filterstring = None
        self.conn = None
        self.cursor = None
        if dbpath != None:
            self.dbpath = dbpath
            if mode == 'sub':
                self.connect_db()
        if filter != None:
            self.filter = filter
        else:
            if filterstring != None:
                self.filterstring = filterstring
            elif filtername != None:
                self.filtername = filtername
            elif filterpath != None:
                self.filterpath = filterpath
            if mode == 'sub':
                self.loadfilter()
        
    def run (self, cmd=None, args=None, dbpath=None, filter=None):
        if args != None:
            self.parse_args(args)
        if cmd != None:
            self.cmd = cmd
        
        if dbpath != None:
            self.dbpath = dbpath
            self.connect_db()
        elif self.dbpath != None and self.cursor == None:
            self.connect_db()
        
        # Loads filter.
        if filter != None:
            self.filter = filter
        elif (self.filtername != None or self.filterpath != None or 
            self.filterstring != None) and self.filter == None:
            self.loadfilter()

        ret = None
        if self.cursor != None and self.filter != None:
            if self.cmd == 'uidpipe':
                ret = self.run_level_based_func(self.getuiditerator)
            elif self.cmd == 'count':
                ret = self.run_level_based_func(self.getcount)
            elif self.cmd == 'rows':
                ret = self.run_level_based_func(self.getrows)
            elif self.cmd == 'pipe':
                ret = self.run_level_based_func(self.getiterator)
        elif self.cursor != None and self.cmd == 'list':
            ret = self.listfilter()
        
        # Saves filter.
        if self.filter != None:
            if self.cmd == 'save' or self.savefiltername != None:
                ret = self.savefilter()
        
        return ret
    
    def run_level_based_func (self, cmd):
        ret = {}
        if self.level != None:
            ret[self.level] = cmd(level=self.level)
        else:
            levels = ['variant', 'gene']
            ret = {}
            for level in levels:
                ret_onelevel = cmd(level=level)
                ret[level] = ret_onelevel
        return ret

    def parse_args (self, args):
        parser = argparse.ArgumentParser()
        parser.add_argument('-d',
            dest='dbpath',
            required=True,
            help='Path of a result database file (.sqlite)')
        parser.add_argument('-f',
            dest='filterpath',
            help='Path of a filtering criteria file')
        parser.add_argument('-F',
            dest='filtername',
            help='Name of the filter to apply (saved in the database)')
        parser.add_argument('--filterstring',
            dest='filterstring',
            default=None,
            help='Filter in JSON')
        parser.add_argument('-l',
            dest='level',
            default=None,
            choices=['variant', 'gene'],
            help='Analysis level to filter')
        parser.add_argument('-s',
            dest='savefiltername',
            help='Name to save the filter as (in the database)')
        if self.mode == 'main':
            parser.add_argument('command', 
            choices=['uidpipe', 'count', 'rows', 'pipe', 'save', 'list'],
            help='Command')
        
        parsed_args = parser.parse_args(args)
        self.dbpath = parsed_args.dbpath
        self.filterpath = parsed_args.filterpath
        self.level = parsed_args.level
        if self.mode == 'main':
            self.cmd = parsed_args.command
        self.savefiltername = parsed_args.savefiltername
        self.filtername = parsed_args.filtername
        self.filterstring = parsed_args.filterstring
        
    def connect_db (self, dbpath=None):
        if dbpath != None:
            self.dbpath = dbpath
        self.conn = sqlite3.connect(self.dbpath)
        self.cursor = self.conn.cursor()
        self.conn.create_function('regexp', 2, regexp)
    
    def close_db (self):
        self.cursor.close()
        self.conn.close()
    
    def create_filtertable (self):
        if self.cursor == None:
            return
        sql = 'create table ' + self.filtertable + ' (name text, criteria text)'
        self.cursor.execute(sql)
        self.conn.commit()
    
    def filtertable_exists (self):
        sql = 'select name from sqlite_master where ' +\
            'type="table" and name="' + self.filtertable + '"'
        self.cursor.execute(sql)
        row = self.cursor.fetchone()
        if row == None:
            ret = False
        else:
            ret = True
        return ret
    
    def loadfilter (self, filterpath=None, filtername=None, 
            filterstring=None, filter=None):
        if filterpath != None:
            self.filterpath = filterpath
        if filtername != None:
            self.filtername = filtername
        if filterstring != None:
            self.filterstring = filterstring
        if filter != None:
            self.filter = filter
        
        if self.filterstring != None:
            self.filterstring = self.filterstring.replace("'", '"')
            self.filter = json.loads(self.filterstring)
        elif self.filtername != None and self.filtertable_exists():
            self.cursor.execute('select criteria from ' + self.filtertable +
                ' where name="' + self.filtername + '"')
            criteria = self.cursor.fetchone()
            if criteria != None:
                self.filter = json.loads(criteria[0])
        elif self.filterpath != None and os.path.exists(self.filterpath):
            with open(self.filterpath) as f:
                ftype = self.filterpath.split('.')[-1]
                if ftype in ['yml','yaml']:
                    self.filter = yaml.load(f)
                elif ftype in ['json']:
                    self.filter = json.load(f)
            
    def delete_filtered_uid_table (self):
        self.cursor.execute('pragma synchronous=0')
        q = 'drop table if exists variant_filtered'
        self.cursor.execute(q)
        q = 'drop table if exists gene_filtered'
        self.cursor.execute(q)
        self.conn.commit()
        self.cursor.execute('pragma synchronous=2')

    def getwhere (self, level):
        keycol = get_keycol(level)
        if self.filter == None:
            q = 'select {} from {}'.format(keycol, level)
        else:
            if level not in self.filter:
                q = 'select {} from {}'.format(keycol, level)
            else:
                criteria = self.filter[level]
                main_group = FilterGroup(criteria, level)
                q = main_group.get_sql()
        return q

    def getvariantcount (self):
        return self.getcount('variant')

    def getgenecount (self):
        return self.getcount('gene')

    def getcount (self, level='variant'):
        level = 'variant'
        self.make_filtered_uid_table()
        ftable = level + '_filtered'
        q = 'select count(*) from ' + ftable
        self.cursor.execute(q)
        n = self.cursor.fetchone()[0]
        if self.stdout == True:
            print('#' + level)
            print(str(n))
        return n

    def getvariantrows (self):
        return self.getrows('variant')

    def getgenerows (self):
        return self.getrows('gene')

    def getrows (self, level='variant'):
        (sample_needed, tag_needed, include_where, exclude_where) = self.getwhere(level)
        q = 'select *  from ' + level + include_where + ' except select * from ' + level + exclude_where
        self.cursor.execute(q)
        ret = [list(v) for v in self.cursor.fetchall()]
        if self.stdout == True:
            print('#' + level)
            for row in ret:
                print('\t'.join([str(v) for v in row]))
        
        return ret

    def get_gene_row (self, hugo):
        q = 'select * from gene where base__hugo=?'
        self.cursor.execute(q, [hugo])
        row = self.cursor.fetchone()
        return row

    def getvariantiterator (self):
        return self.getiterator('variant')
    
    def getgeneiterator (self):
        return self.getiterator('gene')
    
    def getiterator (self, level='variant'):
        (sample_needed, tag_needed, include_where, exclude_where) = self.getwhere(level)
        sql = 'select *  from ' + level + include_where + ' except select * from ' + level + exclude_where
        self.cursor.execute(sql)
        it = self.cursor.fetchall()
        return it

    def get_filtered_iterator (self, level='variant'):
        kcol = get_keycol(level)
        if level == 'variant':
            ftable = 'variant_filtered'
        elif level == 'gene':
            ftable = 'gene_filtered'
        elif level == 'sample':
            ftable = 'variant_filtered'
        elif level == 'mapping':
            ftable = 'variant_filtered'
        table = level
        if level in ['variant', 'gene', 'sample', 'mapping']:
            sql = 'select t.* from ' + table + ' as t inner join ' + ftable +\
                ' as f on t.' + kcol + '=f.' + kcol
        self.cursor.execute(sql)
        it = self.cursor.fetchall()
        return it

    def make_filtered_uid_table (self):
        self.cursor.execute('pragma synchronous=0')
        level = 'variant'
        vtable = level
        vftable = level + '_filtered'
        q = 'drop table if exists ' + vftable
        self.cursor.execute(q)
        q = 'create table {} as select * from '.format(vftable) 
        subq = self.getwhere(level)
        q = q + '(' + subq + ')'
        print(q)
        self.cursor.execute(q)
        self.cursor.execute('pragma synchronous=2')

    def make_filtered_hugo_table (self):
        self.cursor.execute('pragma synchronous=0')
        level = 'gene'
        vtable = 'variant'
        keycol = get_keycol(level)
        vftable = vtable + '_filtered'
        gftable = level + '_filtered'
        q = 'drop table if exists ' + gftable
        self.cursor.execute(q)
        q = 'create table gene_filtered as select distinct v.base__hugo from variant as v inner join variant_filtered as vf on vf.base__uid=v.base__uid where v.base__hugo is not null'
        self.cursor.execute(q)
        self.cursor.execute('pragma synchronous=2')
        
    def savefilter (self, name=None):
        if self.cursor == None or self.filter == None:
            return
        
        if name == None:
            if self.savefiltername != None:
                name = self.savefiltername
            else:
                name = 'default'
        
        # Creates filter save table if not exists.
        self.cursor.execute('select name from sqlite_master where ' +
            'type="table" and name="' + self.filtertable + '"')
        ret = self.cursor.fetchone()
        if ret == None:
            self.cursor.execute('create table ' + self.filtertable +
                ' (name text unique, criteria text)')
        
        # Saves the filter.
        filterstr = json.dumps(self.filter)
        sql = 'insert or replace into ' + self.filtertable +\
            ' values ("' + name + '", \'' + filterstr + '\')'
        self.cursor.execute(sql)
        self.conn.commit()
    
    def listfilter (self, name=None):
        if name == None:
            if self.savefiltername != None:
                name = self.savefiltername
            else:
                name = 'default'
        
        # Creates filter save table if not exists.
        self.cursor.execute('select name from sqlite_master where ' +
            'type="table" and name="' + self.filtertable + '"')
        ret = self.cursor.fetchone()
        if ret == None:
            self.cursor.execute('create table ' + self.filtertable +
                ' (name text, criteria text)')
        
        sql = 'select name, criteria from ' + self.filtertable
        self.cursor.execute(sql)
        ret = {}
        for row in self.cursor.fetchall():
            name = row[0]
            criteria = json.loads(row[1])
            ret[name] = criteria
            if self.stdout:
                print('#' + name)
                for level in criteria:
                    print('    ' + level + ':')
                    for column in criteria[level]:
                        print('        ' + column + ': ' + 
                              criteria[level][column])
            
        return ret
    
    def addvariantfilter (self, column, condition):
        self.addfilter(column, condition, 'variant')
        
    def addgenefilter (self, column, condition):
        self.addcriterion(column, condition, 'gene')
        
    def addfilter (self, column, condition, level='variant'):
        if self.filter == None:
            self.filter = {}
        if level not in self.filter:
            self.filter[level] = {}
        self.filter[level][column] = condition
    
    def removevariantfilter (self, column):
        self.removefilter(column, 'variant')
        
    def removegenefilter (self, column):
        self.removefilter(column, 'gene')
        
    def removefilter (self, column, level='variant'):
        if self.filter == None:
            return
        if level in self.filter and column in self.filter[level]:
            del self.filter[level][column]
    
    def table_exists (self, table):
        sql = 'select name from sqlite_master where type="table" and ' +\
            'name="' + table + '"'
        self.cursor.execute(sql)
        if self.cursor.fetchone() == None:
            return False
        else:
            return True
    
    def get_variant_iterator_filtered_uids_cols (self, cols):
        q = 'select ' + ','.join(cols) + ' from variant as v ' +\
            'inner join variant_filtered as f on v.base__uid=f.base__uid'
        self.cursor.execute(q)
        for row in self.cursor.fetchall():
            d = {}
            for i in range(len(row)):
                d[cols[i].split('__')[1]] = row[i]
            yield d
    
    def get_result_levels (self):
        q = 'select name from sqlite_master where type="table" and ' +\
            'name like "%_header"'
        self.cursor.execute(q)
        table_names = []
        for row in self.cursor.fetchall():
            table_names.append(row[0].replace('_header', ''))
        return table_names

def regexp (y, x, search=re.search):
    if x is None:
        return 0
    return 1 if search(y, x) else 0

def main ():
    cv = CravatFilter(mode='main')
    cv.run(args=sys.argv[1:])

def test ():
    cf = CravatFilter(dbpath='d:\\git\\cravat-newarch\\tmp\\job\\in1000.sqlite',
                      mode='sub',
                      filterstring='{"variant": {"thousandgenomes__af": ">0.1"}}')
    print(cf.getcount())
    cf.make_filtered_uid_table()
    for row in cf.get_filtered_iterator(level='variant'):
        print(row)
    for row in cf.get_filtered_iterator(level='gene'):
        print(row)

if __name__ == '__main__':
    #main()
    test()
