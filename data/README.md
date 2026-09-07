# The Dominick's orange-juice panel

Weekly store-level scanner data on refrigerated orange juice at Dominick's Finer Foods, a
Chicago supermarket chain: 83 stores, 11 brands, 121 weeks (weeks 40 to 160 of the chain's
calendar), 106,139 brand-store-week observations, 96 percent of the complete product. It is the
data set of the empirical application in Harrison, Canavire Bacarreza, Jacho-Chávez and
Rios-Avila (2026), and a genuinely irregular three-way panel: every store carries every brand at
some point and every brand is observed in every week, but 394 of the 10,043 store-weeks are
missing and, within a realized store-week, not every brand is present.

## Getting the data

The files are not shipped with this repository. They are fetched from CRAN and written here by

```stata
cd examples
do 00_get_data.do
```

which requires R with `Rscript` on the path (the script installs the package `bayesm` if it is
missing). `examples/02_orange_juice.do` runs it on its own the first time. Afterwards this folder
holds

| File | Contents |
|---|---|
| `orangeJuice.dta` | one row per (store, brand, week): `logmove` (log of units sold), `price1` to `price11` (shelf price of each brand in the store-week, dollars per ounce), `price` and `lprice` (the row's own brand's price and its log), `deal` (in-store coupon, 0/1), `feat` (share of the week on feature advertisement), `profit` |
| `orangeJuice_storedemo.dta` | the eleven demographic characteristics of each store's trading area, key `STORE` |

Brands: 1 Tropicana Premium 64 oz, 2 Tropicana Premium 96 oz, 3 Florida's Natural 64 oz
(premium tier); 4 Tropicana 64 oz, 5 Minute Maid 64 oz, 6 Minute Maid 96 oz, 7 Citrus Hill
64 oz, 8 Tree Fresh 64 oz, 9 Florida Gold 64 oz (national tier); 10 Dominicks 64 oz,
11 Dominicks 128 oz (store brand).

## Source

The data are the object `orangeJuice` of the R package **bayesm** (Rossi, version 3.1-7),
released under the GNU General Public License, version 2 or later. Their origin is Montgomery,
A. L. (1997), Creating micro-marketing pricing strategies using supermarket scanner data,
*Marketing Science* 16(4): 315–337, who assembled the panel from the Dominick's Finer Foods
database of the **Kilts Center for Marketing, University of Chicago Booth School of Business**,
whose terms ask that the Center be acknowledged in publications using the data. The same data
are used in Chapter 5 of Rossi, Allenby and McCulloch (2005), *Bayesian Statistics and
Marketing*, Wiley.

`00_get_data.do` drops the constant column, adds `price` and `lprice`, and attaches variable
and value labels; the values are unchanged. If you use the data, please cite Montgomery (1997)
as the source, the `bayesm` package as the distribution, and acknowledge the Kilts Center for
Marketing.
