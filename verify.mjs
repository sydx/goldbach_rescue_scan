// Exact Eratosthenes sieve and exhaustive finite checks. No dependencies.
// Run: node verify.mjs > verification.json
const transitionLimit = 10000;
const sixStepLimit = 1000000;
const rescueCrosscheckLimit = 1000000;
const bound = 2 * sixStepLimit + 10;
const prime = new Uint8Array(bound + 1).fill(1);
prime[0] = prime[1] = 0;
for (let d = 2; d * d <= bound; ++d)
  if (prime[d]) for (let x = d*d; x <= bound; x += d) prime[x] = 0;
const odd = [];
for (let p = 3; p <= bound; p += 2) if (prime[p]) odd.push(p);
const isPrime = x => x >= 2 && x <= bound && prime[x] === 1;
function representations(n) {
  const out = [];
  for (const p of odd) {
    if (p > n) break;
    if (isPrime(2*n-p)) out.push([p,2*n-p]);
  }
  return out;
}
let extinctionCount = 0;
const trFailures = [], rightExchangeFailures = [], rightEndpointFailures = [];
const nonModularExtinctions = [];
for (let n = 3; n <= transitionLimit; ++n) {
  const pairs = representations(n);
  const extinct = pairs.length > 0 && pairs.every(([p,q]) =>
    !isPrime(p+2) && !isPrime(q+2));
  if (!extinct) continue;
  ++extinctionCount;
  if (n % 3 !== 1) nonModularExtinctions.push(n);
  if (!pairs.some(([p,q]) =>
    (isPrime(p+4) && isPrime(q-2)) ||
    (isPrime(p-2) && isPrime(q+4)))) trFailures.push(n);
  // q-2 and q must both lie strictly right of the NEW pivot.
  if (!pairs.some(([p,q]) => q-2 > n+1 &&
    isPrime(q-2) && isPrime(p+4))) rightExchangeFailures.push(n);
  const newer = representations(n+1);
  if (!newer.some(([u,v]) => v > n+1 &&
    (isPrime(v+2) || (v-2 > n+1 && isPrime(v-2)))))
    rightEndpointFailures.push(n);
}
const sixBases = odd.filter(p => isPrime(p+6));
let largestLeastP = 0, extremalWitness = null;
const sixFailures = [];
for (let n = 7; n <= sixStepLimit; ++n) {
  const M = 2*n-6;
  let witness = null;
  // Full ordered search: the endpoint with a forward gap six may be larger.
  for (const p of sixBases) {
    if (p > M-3) break;
    if (isPrime(M-p)) { witness = p; break; }
  }
  if (witness === null) sixFailures.push(n);
  else if (witness > largestLeastP) {
    largestLeastP = witness;
    extremalWitness = {n,p:witness,q:M-witness,pPlus6:witness+6};
  }
}
const at277 = representations(277);
const rightCandidates277 = at277.filter(([p,q]) =>
  q-2 > 278 && isPrime(q-2)).map(([p,q]) =>
  ({old:[p,q],candidate:[p+4,q-2],leftIsPrime:isPrime(p+4)}));
// A second, dense-byte sieve implementation of the rescue search, scanning
// cousin lower members rather than the C++ program's twin-first shortcut.
const twins = odd.filter(t => isPrime(t+2));
const cousins = odd.filter(c => isPrime(c+4));
const rescueFailures = [], otherExtinctions = [];
let residueTargets = 0, rescueCrosscheckMaxC = 0;
function hasRescue(n) {
  const M = 2*n-2;
  for (const c of cousins) {
    if (c > M-3) break;
    const t = M-c;
    if (isPrime(t) && isPrime(t+2)) {
      rescueCrosscheckMaxC = Math.max(rescueCrosscheckMaxC,c);
      return true;
    }
  }
  return false;
}
function hasDirectTwin(n) {
  const N = 2*n;
  for (const t of twins) {
    if (t > N-3) break;
    if (isPrime(N-t)) return true;
  }
  return false;
}
for (let n=4; n<=rescueCrosscheckLimit; ++n) {
  if (n%3===1) {
    if (isPrime(2*n-3)) continue;
    ++residueTargets;
    if (!hasRescue(n)) rescueFailures.push(n);
  } else {
    if (hasDirectTwin(n)) continue;
    if (representations(n).length) {
      otherExtinctions.push(n);
      if (!hasRescue(n)) rescueFailures.push(n);
    }
  }
}
const results = {
  method: 'Exact Eratosthenes sieve; exhaustive integer ranges',
  sieveBound: bound,
  transitions: {range:[3,transitionLimit],extinctionCount,trFailures,
    rightExchangeFailures,rightEndpointFailures,nonModularExtinctions},
  sixStep: {range:[7,sixStepLimit],tested:sixStepLimit-6,
    failures:sixFailures,largestLeastP,extremalWitness},
  independentRescueCrosscheck: {range:[4,rescueCrosscheckLimit],
    residueTargets,rescueFailures,otherExtinctions,
    greatestLeastCousinBase:rescueCrosscheckMaxC},
  firstDirectionalCounterexample: {n:277,old:at277,
    rightCandidates:rightCandidates277,new:representations(278)}
};
if (extinctionCount !== 2197 || trFailures.length ||
    JSON.stringify(rightExchangeFailures) !== '[277,676,1012,1342]' ||
    rightEndpointFailures.length || sixFailures.length || largestLeastP !== 2281 ||
    residueTargets !== 258814 || rescueFailures.length || otherExtinctions.length)
  throw new Error('Recomputed results differ from the manuscript; inspect results.');
console.log(JSON.stringify(results,null,2));
