#!/bin/bash
for i in $(ls |grep sys2 | grep xvg ) 
do case "$i" in 
	*-1.xvg) mv $i $( echo $i | sed -e 's/-1.xvg/-bond.xvg/g' ) ;;
	*-2.xvg) mv $i $( echo $i | sed -e 's/-2.xvg/-angle/g') ;;
	*-3.xvg) mv $i $( echo $i | sed -e 's/-3.xvg/-propdih.xvg/g' ) ;;
	*-4.xvg) mv $i $( echo $i | sed -e 's/-4.xvg/-impdih.xvg/g') ;;
	*-5.xvg) mv $i $( echo $i | sed -e 's/-5.xvg/-lj14.xvg/g') ;;
	*-6.xvg) mv $i $( echo $i | sed -e 's/-6.xvg/-coulomb14.xvg/g') ;;
	*-7.xvg) mv $i $( echo $i | sed -e 's/-7.xvg/-ljsr.xvg/g' ) ;;
	*-8.xvg) mv $i $( echo $i | sed -e 's/-8.xvg/-dispcorr.xvg/g') ;;
	*-9.xvg) mv $i $( echo $i | sed -e 's/-9.xvg/-coulombsr.xvg/g' ) ;;
	*-10.xvg) mv $i $( echo $i | sed -e 's/-10.xvg/-coulrecip.xvg/g') ;;
	*-11.xvg) mv $i $( echo $i | sed -e 's/-11.xvg/-potential.xvg/g') ;;
esac
done
