import { nucleusToElement } from './smiles-to-hose.js';

describe('nucleusToElement', () => {
  test('13C -> C', () => {
    expect(nucleusToElement('13C')).toBe('C');
  });

  test('1H -> H', () => {
    expect(nucleusToElement('1H')).toBe('H');
  });

  test('15N -> N', () => {
    expect(nucleusToElement('15N')).toBe('N');
  });

  test('31P -> P', () => {
    expect(nucleusToElement('31P')).toBe('P');
  });

  test('19F -> F', () => {
    expect(nucleusToElement('19F')).toBe('F');
  });

  test('11B -> B', () => {
    expect(nucleusToElement('11B')).toBe('B');
  });

  test('29Si -> Si', () => {
    expect(nucleusToElement('29Si')).toBe('Si');
  });
});
