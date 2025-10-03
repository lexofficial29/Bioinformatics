const fs = require('fs');
const readline = require('readline');

const FILENAME = process.argv[2];

if (!FILENAME) {
  console.error("Please provide a FASTA filename as the first argument.");
  process.exit(1);
}

let alphabet = {};
let total = 0;

const rl = readline.createInterface({
  input: fs.createReadStream(FILENAME),
  crlfDelay: Infinity
});

rl.on('line', (line) => {
  const trimmed = line.trim();
  if (trimmed.startsWith('>')) return;

  for (let char of trimmed.toUpperCase()) {
    if (!'ACGTN'.includes(char)) continue;

    if (!alphabet[char]) {
      alphabet[char] = 1;
    } else {
      alphabet[char] += 1;
    }
    total += 1;
  }
});

rl.on('close', () => {
  let frequency = {};
  for (let key in alphabet) {
    frequency[key] = ((alphabet[key] / total) * 100).toFixed(2) + "%";
  }

  console.log("Frequencies:", frequency);
});
