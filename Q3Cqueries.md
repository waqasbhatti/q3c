**The '''cone search''' (the query of all objects within the circular region of the sky).**

To query all objects lying within radius0=0.1 deg from ra0 = 11deg,
dec0 = 12deg in table mytable you should do:
```
   SELECT * FROM mytable WHERE q3c_radial_query(ra, dec, 11, 12, 0.1);
```
IMPORTANT limitation is that those queries are only efficient is the area is not too large. E.g if the diameter is ~ 40-50 degrees or so those queries will be inefficient

**The '''polygonal query''' (the query of the objects which lie inside the region
bounded by the polygon on the sphere)**

To query the objects in the polygon ((0,0),(2,0),(2,1),(0,1)) )
(this is the spherical polygon with following vertices:
(ra=0, dec=0) ; (ra=2, dec=0); (ra=2, dec=1); (ra=0, dec=1)):
```
  SELECT * FROM mytable WHERE
        q3c_poly_query(ra, dec, '{0, 0, 2, 0, 2, 1, 0, 1}');
```
**The '''ellipse query''' (the query of objects which lie inside an ellipse, with
ra=10, dec=20, major axis=1deg, axis ratio=0.5 and positional angle=10deg)
```
   SELECT * FROM mytable WHERE
    q3c_ellipse_query(ra, dec, 10, 20, 1, 0.5 ,10);
```**


**The positional''' cross-match''' of the tables:**

Lets Assume that we have a huge table "table2" with ra and dec columns and
an already created Q3C index and the smaller table "table1" with ra and dec
columns.

Now, if we want to cross-match the tables "table1" and "table2" by position
with the crossmatch radius of 0.0002 degrees, we should do:
```
   SELECT * FROM table1 AS a, table2 AS b WHERE
        q3c_join(a.ra, a.dec, b.ra, b.dec, 0.0002);
```
If every object in the table1 have his own error circle (like with X-ray data
for example) (suppose that the radius of that circle is the column "err"),
then you can run the query:
```
   SELECT * FROM table1 AS a, table2 AS b WHERE
        q3c_join(a.ra, a.dec, b.ra, b.dec, a.err);
```
'''Important:''' If additionally to the crossmatch you want to filter some rows by some condition, the query will become more complicated, e.g.:
```
   WITH a as (
             SELECT * from table1 where column_from_table1>10
             )
   SELECT * FROM a, table2 AS b WHERE
         q3c_join(a.ra, a.dec, b.ra, b.dec, 0.0002);
```
This will ensure that the crossmatch will execute efficiently;

If you want to also filter our some rows using the columns from table2, the query becomes more complex

```
  WITH a AS (
             SELECT * from table1 where column_from_table1>10
             ),
  ab AS (
             SELECT * from a, table2 AS b WHERE
                   q3c_join(a.ra, a.dec, b.ra, b.dec, 0.0002);
             )
   SELECT * FROM ab WHERE column_from_table2<44
```

**The positional cross-match of the tables with the ellipse error-area:
(for example if you want to find all the objects from one catalogue which lies
inside the elliptical bodies of the galaxies from the second catalogue)**

It is possible to do the join when the error area of each record of the
catalogue is an ellipse. Then you can do the query like this
```
 SELECT * FROM table1 AS a, table2 AS b WHERE
        q3c_ellipse_join(a.ra, a.dec, b.ra, b.dec, a.maj_ax
        a.axis_ratio, a.PA);
```
where axis\_ratio is the column with axis ratio of the ellipses and PA is the
column with the positional angles of them, and maj\_ax is the column with major
axises of those ellipses.


**If you want to find sources in the first table and doesn't have a match in the second table, the way to do it is:
```
  SELECT * FROM table1 as a where (
              select count(*) from table 2as b where q3c_join(a.ra,a.dec,b.ra,b.dec, 0.0002) ) = 0;
```**

**'''Nearest neighbor''' (works in PG 9.3+)**

This query selects the only nearest neighbor for each row in your table.
If there are no neighbor, the columns are filled with nulls.
```
 SELECT  t.*, ss.* FROM mytable AS t,
        LATERAL (
                SELECT s.* 
                     FROM 
                         sdssdr9.phototag AS s
                     WHERE
                         q3c_join(t.ra, t.dec, s.ra, s.dec, 1./3600)
                     ORDER BY
                         q3c_dist(t.ra,t.dec,s.ra,s.dec)
                     ASC LIMIT 1
                ) as ss;
```
The idea is very simple for every row of your table mytable LATERAL() executes the "subquery" orders them by distance and limit by 1.

**'''Nearest neighbor 2''' (works with PG 9.0+)**

This query selects the only nearest neighbor for each row in your table.
If there are no neighbor, the columns are filled with nulls. This query requires presence of some object id column with the index on the table.
```
 WITH x AS (
       SELECT *, ( SELECT objid FROM sdssdr9.phototag AS p WHERE q3c_join(m.ra, m.dec, p.ra, p.dec, 1./3600)
                   ORDER BY q3c_dist(m.ra, m.dec, p.ra, p.dec) ASC LIMIT 1) AS match_objid  FROM mytable AS m 
           )
 SELECT * FROM x, sdssdr9.phototag AS s WHERE x.match_objid=s.objid;
```

**Angular correlation function**

```

SELECT width_bucket(q3c_dist(a.ra, a.dec, b.ra, b.dec), 0, 1, 100) AS bin,
       count(*) AS cnt 
       FROM twoxmm.main AS a,
            twoxmm.main AS b 
       WHERE q3c_join(a.ra, a.dec, b.ra, b.dec, 1) AND 
             a.srcid != b.srcid 
       GROUP BY width_bucket(q3c_dist(a.ra, a.dec, b.ra, b.dec), 0, 1, 100);

```