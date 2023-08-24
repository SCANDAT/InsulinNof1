%let localpath=%substr(%sysget(SAS_EXECFILEPATH),1,%eval(%length(%sysget(SAS_EXECFILEPATH))-%length(%sysget(SAS_EXECFILENAME))-1));;
libname nof1 "&localpath\..\Data";
*Import CGM data;
data allcgm;
infile "&localpath\..\Data\CGM.csv" delimiter = ',' MISSOVER DSD lrecl=13106 firstobs=2 ;
   informat Zulu_Time anydtdtm40. ;
   informat Local_Time anydtdtm40. ;
   informat Value best32. ;
   informat Units $6. ;
   format Zulu_Time datetime. ;
   format Local_Time datetime. ;
   format Value best12. ;
   format Units $6. ;
input
            Zulu_Time
            Local_Time
            Value
            Units  $
;

run;
*Summarize at 5-minute intervals;
proc sql;
create table allcgm2 as
select distinct
  300*floor(local_time/300) as local_time format=datetime.,
  mean(value) as value
from allcgm
group by 300*floor(local_time/300);
quit;
*Sort and make sure no duplicates;
proc sort data=allcgm2 nodupkey force;
by Local_Time value;
run;

*Import the protocol;
data protocol;
infile "&localpath\..\Data\Protocol.csv" delimiter = ',' MISSOVER DSD lrecl=13106 firstobs=2 ;
informat Block best32. ;
informat insulintype $1. ;
informat Startdatum yymmdd10. ;
informat Starttid time20.3 ;
informat Slutdatum yymmdd10. ;
informat Sluttid time20.3 ;
format Block best12. ;
format insulintype $1. ;
format Startdatum yymmdd10. ;
format Starttid time20.3 ;
format Slutdatum yymmdd10. ;
format Sluttid time20.3 ;
input
         Block
         insulintype  $
         Startdatum
         Starttid
         Slutdatum
         Sluttid
;
format startdt enddt datetime.;
startdt=dhms(startdatum,0,0,starttid);
enddt=dhms(Slutdatum,0,0,Sluttid);
time=(enddt-startdt)/3600;
run;
*Check that time makes sense;
proc means data=protocol;
class insulintype;
ways 0 1;
var time;
run;
*Combine to create analysis dataset;
*In main analysis, no washout;
proc sql;
create table studydataset as
select distinct
  a.*,
  b.local_time,
  b.value
from protocol a inner join allcgm2 b
  on a.startdt lt b.local_time lt a.enddt
order by b.local_time;
quit;
proc means data=studydataset;
class insulintype;
var value;
run;

*Store for further analysis with R and define high vs low values;
data nof1.studydataset;
set studydataset;
above=(value gt 10);
person=1;
below=(value lt 4);
run;
*Summarize basics;
proc freq data=nof1.studydataset;
tables insulintype*(above below) / nocol nopercent;
run;

proc sort data=nof1.studydataset;
by local_time;
run;
data sd3;
set nof1.studydataset;
nr=_n_;
run;
*This shit takes forever to run, so don't try it;
proc hpmixed data=sd3;
class insulintype nr person;
model value=insulintype / s;
repeated nr / type=ar(1) subject=person;
run;

*Here we analyze risk of high vs low values comparing the two treatments;
proc genmod data=nof1.studydataset;
class insulintype person;
model above(event='1')=insulintype / dist=binomial link=logit ;
repeated subject=person /type=ar(1) modelse corrb covb;
run;


proc genmod data=nof1.studydataset;
class insulintype person;
model below(event='1')=insulintype / dist=binomial link=logit;
repeated subject=person /type=ar(1) modelse corrb covb;
run;

*Explore a 2-hour washout (no big change);

proc sql;
create table sdwashout as
select distinct
  a.*,
  b.local_time,
  b.value
from protocol a inner join allcgm2 b
  on a.startdt+7200 lt b.local_time lt a.enddt
order by b.local_time;
quit;


data nof1.studydataset_wo;
set sdwashout;
above=(value gt 10);
person=1;*"G";
below=(value lt 4);
run;

proc freq data=nof1.studydataset_wo;
tables insulintype*(above below) / nocol nopercent;
run;

proc means data=nof1.studydataset_wo mean median mode p99 max;
class insulintype;
var value;
run;


proc genmod data=nof1.studydataset_wo;
class insulintype person;
model above(event='1')=insulintype / dist=binomial link=logit ;
repeated subject=person /type=ar(1) modelse corrb covb;
lsmestimate insulintype 1 -1 /exp cl;
run;


proc genmod data=nof1.studydataset_wo;
class insulintype person;
model below(event='1')=insulintype / dist=binomial link=logit;
repeated subject=person /type=ar(1) modelse corrb covb ;
lsmestimate insulintype 1 -1 /exp cl;
run;
