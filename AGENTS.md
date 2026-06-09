# Conventions for continuing this project

This file is the working agreement for any person or model continuing the work. The sibling project
established these through repeated correction; they are not optional polish.

## Tone and writing

- Matter-of-fact and precise. State what was done and what was found. No hype, no slogans, no
  congratulatory or motivational framing.
- Do not oversell a result and then walk it back. If a claim is uncertain, say so up front.
- When something is wrong or a mistake was made, say so plainly in one line and move on. Do not frame
  a failure or a self-correction as if the error were valuable.
- Be brief. Prefer concrete numbers and file references over abstraction.
- No em-dashes. No first or second person in the public docs.
- Factual figure titles. No figure or claim that the code and results do not support.

## Honesty and verification

- Every quantitative claim must be backed by code and a result that can be regenerated. A number that
  appears only in prose is a defect.
- Beating chance, or beating a baseline that cannot do the task, is not a result. State the real
  comparator.
- Verify the method, the code, and the write-up before reporting. Re-running and finding an error is
  expected, not noteworthy.

## Calibration and uncertainty in claims

- Frame findings as hypotheses with the evidence for and against, not as settled conclusions. State
  the leading explanation, why it is credible, and the reasons it might be wrong.
- Do not claim full understanding of why a result occurred. Give the most likely cause and the
  alternatives that have not been ruled out.
- Do not assert an exact or complete solution. Report the actual numbers and what they do and do not
  establish.
- Avoid self-certifying words ("honest", "honestly", "to be clear", "the real answer"). They add no
  information and imply the rest needs certifying.
- Stay calibrated to the evidence. Do not alternate between overselling and underselling. A result is
  usually partial; saying how partial it is, is the point.
- When proposing work, name the problem first, then the candidate approaches, then why each is or is
  not promising. Prioritize among them rather than presenting one option as the answer.

## Privacy (this repo is public-facing)

- Citing published papers by author and year is fine and expected (see `docs/literature.md`).
- Do not state or imply collaboration with, mentorship by, or endorsement from any individual.
- Do not name any individual's employer, institution, school, or role.
- Keep raw data, private notes, and personal context out of the public tree.

## Working style

- Default to exploring more than one option and prioritizing among them, rather than forcing a single
  blind choice. Say which is prioritized and why.
- The project will evolve from one concrete target outward as more ideas are tested. Keep the goals,
  research questions, and roadmap current as that happens (`docs/update_map.md`).
- Transcribe published equations and parameters from the source, never from memory.

## Git

- Commit and push only when asked. If on the default branch, branch first for non-trivial work.
- Do not add a co-author trailer to commits.
- Do not skip hooks or bypass signing.
