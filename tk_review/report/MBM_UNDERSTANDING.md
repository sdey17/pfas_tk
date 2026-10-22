# What this subject is actually about

*A plain-prose companion to `MBM_PRIMER.md`. No equations to speak of. The
primer has the machinery; this has the ideas. Read this first.*

---

## 1. A number that isn't what it seems

A toxicologist puts liver cells in a plastic dish, adds a chemical at ten
micromolar, waits a day, and writes down what happened. The dose, as reported
in the paper, is ten micromolar.

The cells never saw ten micromolar.

From the moment the chemical enters the dish it starts going places. Some of it
sticks to the plastic walls. Some dissolves into the fatty droplets that come
with the serum in the growth medium. A great deal of it — for the chemicals we
care about here, almost all of it — latches onto albumin, the carrier protein
that serum is full of. Only what is left over, still dissolved in the water,
is free to drift into a cell. By the time things settle down, the free
concentration can be a thousandth of what was weighed out.

This matters enormously, because the reason anyone runs the experiment is to
say something about people. The chain of reasoning goes: *this concentration
harmed cells in a dish, therefore this blood level would be concerning in a
person.* But a dish is not a bloodstream. A dish might be ten percent serum or
half a percent; blood is blood. If you compare the number written on the dish
with a number measured in blood, you are comparing two quantities that were
never the same kind of thing.

The whole field exists to fix that comparison. Everything that follows is
machinery in service of one question: **of everything you put in the well, how
much actually reached the cell?**

---

## 2. The idea: treat the dish as a budget

The fix is almost embarrassingly simple, and that simplicity is worth dwelling
on, because the rest of the subject is just consequences of it.

Think of the well as a set of compartments: the water, the albumin, the fatty
material, the cells themselves, the plastic. You know exactly how much
chemical went in. None of it is destroyed — PFAS in particular are famously
indestructible — and none escapes, because these molecules do not evaporate.
So the amount you put in must still be there, distributed among the
compartments. It is a budget, and budgets balance.

To work out the split you need only one thing per compartment: how strongly
that compartment holds the chemical, compared with plain water. Give a number
to each of those preferences, and the arithmetic of "the parts sum to the
whole" does the rest. Out falls the free concentration in the water, and from
that, everything else.

That is the entire model. One accounting identity and a list of preferences.
When you read that two research groups have published competing mass-balance
models, they are not disagreeing about the accounting. They are disagreeing
about the list.

---

## 3. The quantity that does all the work

The preference number has a name — a partition coefficient — and the easiest
way to understand it is to picture the oldest experiment in the subject. Put
some oil and some water in a jar. Add a pinch of your chemical. Shake it, let
it settle into two layers, and measure the concentration in each. The ratio of
the two is the partition coefficient. It says, in one number, how much this
chemical would rather be in oil than in water.

Two features of that number matter here, and both come back to haunt us.

**It is a ratio, not an amount.** A chemical with a large partition coefficient
for albumin is not necessarily mostly bound to albumin — that also depends on
how much albumin is present. A strong preference for a phase that barely exists
amounts to nothing. This sounds obvious written down, and it is the source of
one of the two large errors we found in this project.

**It is only a constant while the chemical is dilute.** The underlying physics
says the ratio holds when each molecule behaves as if it were alone, unaware of
the others. Crowd enough of them in and they start to interfere — the oil gets
saturated, the binding pockets on the protein fill up — and the "constant"
quietly stops being constant. Nearly every subtlety in this subject is a
version of that assumption failing.

---

## 4. Why PFAS break the mould

These models were built for a particular kind of chemical: neutral, greasy
things like pesticides and industrial solvents, which partition into fat in a
way that tracks hydrophobicity very well. For that family, one measurement
predicts everything. PFAS are not in that family, for three reasons.

**They are surfactants.** A PFAS molecule is a water-repelling fluorinated tail
attached to a water-loving acidic head — the same architecture as soap. Soap
does not behave like oil. It collects at surfaces, forms films, and does things
that a single "greasiness" number cannot anticipate.

**They are charged at body pH.** PFAS are strong acids, so in blood or in
culture medium essentially every molecule carries a negative charge. Charge
changes a molecule's preferences drastically — a charged species is far happier
in water than its neutral twin. The models handle this by applying a penalty to
the neutral prediction, but the size of that penalty is a modelling choice, and
for albumin specifically, one of the models decides the penalty is *zero* —
that the charged form binds albumin exactly as well as the neutral form would.
That single choice is a large part of why albumin ends up dominating these
calculations.

**They bind albumin in a different way entirely.** Most chemicals dissolve into
a protein the way sugar dissolves into tea — diffusely, with no particular
destination. PFAS dock into specific pockets on albumin, the same pockets that
evolved to carry fatty acids. Pockets are countable, and countable things fill
up. One of the two models is built entirely around that fact; the other ignores
it.

---

## 5. Two philosophies

Set side by side, the two published models are a clean illustration of two ways
to build any model at all.

**Armitage's approach is to predict everything from one number.** If you know
how a chemical splits between oil and water, you can estimate how it splits
between albumin and water, between cell membrane and water, between plastic and
water — each by its own calibrated relationship. One input per chemical, and
the whole well falls out. The appeal is enormous: it works for any chemical
anyone has ever measured a greasiness number for, which is tens of thousands of
them, and that is exactly what you need if your job is to screen a national
inventory of industrial chemicals. It was built for breadth.

**Qin's approach is to go to the bench.** Measure how PFAS actually bind
albumin, across a range of concentrations. Find — because this is what you do
find — that the binding curve is not a straight line. At low concentrations
each additional molecule binds as strongly as the last; at high concentrations
the pockets are filling, and the next molecule has only weak, non-specific
sticking available to it. Build that curve into the model so the albumin term
weakens as you load it up. The result is far more faithful to what albumin
does, but it exists only for the handful of chemicals that were measured.
Everything else has to be guessed by interpolation. It was built for depth.

Neither choice is wrong. They are answers to different questions, and the
tension between them runs through everything below.

---

## 6. What is a cell made of?

Here the two models part company completely, and it is worth seeing how far
apart they get.

Armitage builds a cell as a recipe: so much storage fat, so much membrane fat,
so much protein, the rest water, each ingredient contributing its own
preference in proportion to how much of it there is. For a liver-derived cell
line the recipe is overwhelmingly water — about ninety-four percent — with
only a trace of fat and a few percent protein.

But the convention Armitage inherits says that most of a cell's protein does
not sorb chemicals at all; only a small percentage of it behaves like a
chemically interesting phase. Applied to PFAS, that assumption has a dramatic
consequence. For a short-chain PFAS, which has little interest in fat, the
model's cell is essentially a bag of water. For a long-chain one, which loves
fat, the tiny membrane-fat fraction comes to dominate everything, and the model
cell becomes, in effect, a droplet of lipid. The protein contributes a fraction
of one percent throughout. **The convention has defined protein out of the
cell.**

Qin's cell is the opposite. It is protein — over ninety percent of its sorptive
capacity, at every chain length — with lipid playing a supporting role and no
water phase at all. Where Armitage sees a water droplet turning into a fat
droplet, Qin sees a lump of protein.

These are not small differences of emphasis. They are different pictures of
what a cell is. And here is the thing that should make you uneasy: *the two
models' final answers often agree anyway*. They agree because the albumin in
the medium is so dominant that it swamps the cell term, and both models are
really doing an albumin calculation with a cell attached. Agreement between
two models is not corroboration of either when a third thing is doing the work.

---

## 7. Three experiments that do not agree with each other

Before you can say a model is wrong you need data to be wrong against. Three
groups have put PFAS on liver cells and measured how much got inside.

They used different cells — two used a standard liver cancer line, the third
used a more sophisticated line coaxed into behaving like real liver tissue.
They used different amounts of serum, which as we have seen controls how much
albumin is competing with the cells. They dissolved their chemicals in
different solvents at wildly different strengths. They used different plate
formats, so the ratio of plastic surface to liquid differs several-fold. One of
them reported its results as a fraction of what was added; another as a
fraction of what was *recovered*, which is not the same denominator.

And their results disagree. For one chemical at one dose, the third study
reports roughly fourteen times more uptake than the other two. Even the two
studies using the *same* cell line disagree by a factor of three or four on the
short-chain compounds.

That is the honest starting point. A model that disagrees with these data by a
factor of four is disagreeing with a target that is itself uncertain by a
factor of four. Any conclusion has to survive that.

---

## 8. A sequence of disappointments

What happened when we tested the two models is a story worth telling in order,
because the shape of it is more instructive than any individual result.

**Both models over-predict.** They both say more chemical should end up inside
the cells than the experiments find. Not by a little.

**So we looked for the wrong parameter.** The obvious suspect was albumin,
because albumin carries almost the entire budget — well over ninety percent for
the longer chains. If one term holds ninety-five percent of the answer, an
error there is an error in the answer. There is a very good modern measurement
of PFAS–albumin binding, done at realistic concentrations, so we substituted it.

**It made both models worse.** The measurement says albumin binds PFAS *less*
tightly than the models assume. A weaker albumin means less competition for the
chemical, which means more of it free, which means more of it in the cell — and
the models already predicted too much.

**So we tried the cell instead.** There is a measurement of how strongly PFAS
bind cell protein, which is the term Armitage's convention had nearly erased.
Substituting it made both models worse too — because that measurement says cell
protein binds *more* than assumed, which enlarges the cell's appetite, which
again sends more chemical inside.

**Every correction pushed the same way.** This is the pattern worth taking from
the whole exercise: the best available measurement of each parameter moved the
prediction further from the data, every time, because all the measurements
point the same direction. You can be confident a parameter matters, confirm it
is the most influential thing in your model, and still find it is not the thing
that is wrong. **Sensitivity is not error.**

**Then we found a mistake of our own.** Reading the methods sections properly —
which should have happened first — turned up two sentences saying that one of
the three experiments had been run in a serum-*poor* medium, five times leaner
than the growth medium quoted elsewhere in the paper. Our model had been giving
that experiment five times more albumin than it had.

Fixing it made the fit *worse*, and that is the interesting part. The dataset
that had looked like the models' best case was only looking good because of
our error. Remove the error and it joins the others in disagreeing badly. A
model that looks right for the wrong reason is more dangerous than one that
looks wrong, because nobody investigates it.

**Only one thing helps.** The size of a cell. The published figures for how
large a liver cell is span a factor of three, none of the three experiments
measured their own, and adjusting it moves everything by roughly half the
remaining gap. It is the one parameter that rescales every chemical by the same
amount, which is what an unexplained offset needs. It is also the one parameter
nobody bothered to measure, which tells you something about what the field
thinks is interesting.

---

## 9. The disagreement that depends on chain length

One pattern survived all of this and is worth understanding, because it says
something structural rather than accidental.

PFAS come in families that differ by the length of their fluorinated tail. The
longer the tail, the more the molecule loves fat and the more tightly it binds
albumin. Armitage's over-prediction is not constant across that family — it
grows steadily as the chains get longer, reaching about ten-fold for the
longest. Qin's does not grow at all.

The explanation turns out to be elegant. Because albumin dominates the budget,
what controls how much ends up in a cell is not the cell's appetite on its own,
nor albumin's on its own, but the *ratio* between them. Whichever grows faster
as you lengthen the chain wins.

Qin has that ratio roughly constant: his cell and his albumin grow at nearly
the same rate, so the balance between them holds steady and there is no drift.
Armitage has it growing, for two reasons that compound. His albumin
relationship flattens out for very greasy chemicals — a deliberate feature,
calibrated on ordinary drugs, where binding pockets do eventually stop
rewarding extra greasiness. And his cell, being effectively a droplet of fat,
keeps climbing steeply. A fat droplet gets more attractive with chain length
faster than albumin does, so the ratio drifts, and the error grows.

So the two models fail in genuinely different ways. **Qin has the chemistry of
chain length roughly right and the overall level wrong. Armitage has the level
closer and the chain-length behaviour wrong.** They are not two attempts at the
same thing with different accuracy; they are wrong about different things.

---

## 10. The problem underneath all of it

Here is the thing that, once seen, reframes everything above.

How crowded is the albumin in these experiments? Each albumin molecule offers
roughly one good binding site. In the dishes, there is somewhere between a
third of a PFAS molecule and fourteen PFAS molecules for every albumin
molecule. In human blood, there is about one PFAS molecule for every two
thousand albumin molecules.

The experiments are running three to four orders of magnitude more crowded than
the situation they are meant to inform. Their albumin is saturated; the body's
albumin is essentially empty. And as we saw in section 3, a partition
coefficient is only a constant while things are uncrowded. In these dishes it
is not a constant at all — over a single dose series the apparent strength of
albumin binding changes sevenfold.

This is why Qin's model needed its extra machinery: he measured the curve, saw
it bend, and built the bending in. Armitage's model assumes a straight line,
which is defensible for the dilute neutral chemicals it was designed for and is
simply not true here.

But the deeper implication is uncomfortable, and it is the honest conclusion of
this whole line of work. Even a *perfect* model of these dishes would be
describing a chemical regime that does not occur in a person. Better
bookkeeping inside the well gives you a better description of the experiment.
It does not obviously give you a better answer about human exposure, because
the experiment is in the wrong regime to begin with.

---

## 11. The same idea, somewhere it works

It would be easy to conclude from all this that the mass-balance idea is
unsound. It isn't, and there is a clean demonstration.

Take the identical reasoning to blood serum. Serum is a simpler system than a
culture well: water, albumin, and the globulins, with no plastic and no cells.
One group measured, separately and by the same method, how strongly two dozen
PFAS bind to purified albumin, how strongly they bind to purified globulin, and
how strongly they bind to whole serum. That is exactly the test: predict the
third from the first two, with nothing fitted.

It works. The predictions sit on the measurements with no systematic bias and a
typical error of about fifty percent — which, for a quantity spanning a
hundred-fold across the compound set, is good. Letting the protein
concentrations float to improve the fit returns the physiological values and
improves nothing, which is the strongest evidence of all that the structure is
right rather than being propped up by its parameters.

So the accounting is sound. Free concentration times a sum of appetites really
does describe a protein-rich fluid. What fails is the *dish* — a system with
more phases, more uncertain parameters, a plastic term nobody can pin down,
cells whose size nobody measured, and a crowding regime where the central
assumption does not hold.

---

## 12. What the subject teaches

Strip away the chemistry and several lessons generalise well beyond PFAS.

**A model is an accounting identity plus a list of preferences.** The identity
is never in doubt. All the argument is about the list, and about whether each
preference was measured in the regime where it is being used.

**Sensitivity is not error.** Albumin controls almost the entire answer in these
wells, which makes it the obvious suspect and makes it tempting to tune. It is
not where the error is. The two are different questions and conflating them
wastes a great deal of time.

**An error that scales everything cancels; an error that changes a ratio
survives.** The three labs used solvents differing thirty-fold in strength, and
it hardly matters, because a solvent weakens every preference by about the same
factor and the model's answer depends on their ratios. The serum error, which
changed how much albumin there was *relative to* everything else, was worth a
factor of four. Knowing which kind of error you are looking at tells you
immediately whether to care.

**Agreement is not corroboration.** The two models hold opposite views of what a
cell is and frequently agree on the answer, because a third term dominates
both. Two methods agreeing means something only if they could have disagreed.

**The methods section is not optional.** Two sentences about serum content,
buried in a paragraph about cell culture, were worth more than every parameter
substitution we tried. Nobody writes papers about reading the methods section
carefully, and it was the highest-yield thing done in this project.

**Know which regime your measurement came from.** A binding constant measured on
saturated protein and a binding constant measured on empty protein are
different numbers for the same quantity, differing by more than an order of
magnitude. The literature mostly does not record which was which, which makes a
great deal of it unusable without detective work.

---

## 13. When you want the machinery

Everything asserted here is derived, with numbers, in the companion documents:

- `MBM_PRIMER.md` — the equations, the approximations named one by one, the
  derivations, and worked numerical checks. Read it when you want to see *why*
  something above is true rather than take it on trust.
- `MBM_Qin_vs_Armitage.md` — the actual comparison against the three datasets,
  the medium correction, and the chain-length analysis.
- `ALBUMIN_BINDING.md` — albumin binding in its own right, including the serum
  result of section 11.
