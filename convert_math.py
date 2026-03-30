#!/usr/bin/env python3
"""Convert Unicode math symbols to LaTeX in Markdown files.

ONLY handles inline replacements outside fenced code blocks.
Does NOT modify content inside fenced code blocks, inline code, or existing LaTeX.
"""

import re
import sys

GREEK = {
    'ρ': r'\rho', 'α': r'\alpha', 'ν': r'\nu', 'μ': r'\mu',
    'ε': r'\varepsilon', 'ω': r'\omega', 'σ': r'\sigma', 'κ': r'\kappa',
    'γ': r'\gamma', 'δ': r'\delta', 'τ': r'\tau', 'θ': r'\theta',
    'η': r'\eta', 'φ': r'\varphi', 'β': r'\beta', 'λ': r'\lambda',
    'π': r'\pi', 'ζ': r'\zeta',
    'Δ': r'\Delta', 'Ω': r'\Omega', 'Γ': r'\Gamma', 'Λ': r'\Lambda',
    'Σ': r'\Sigma', 'Π': r'\Pi', 'Φ': r'\Phi', 'Ψ': r'\Psi',
    'Θ': r'\Theta', 'Ξ': r'\Xi',
}

OPERATORS = {
    '∇': r'\nabla', '∂': r'\partial', '×': r'\times', '≈': r'\approx',
    '÷': r'\div', '±': r'\pm', '≠': r'\neq', '≤': r'\leq', '≥': r'\geq',
}

SUPER = {
    '⁰':'0','¹':'1','²':'2','³':'3','⁴':'4',
    '⁵':'5','⁶':'6','⁷':'7','⁸':'8','⁹':'9','⁻':'-','⁺':'+',
}

SUB = {
    '₀':'0','₁':'1','₂':'2','₃':'3','₄':'4',
    '₅':'5','₆':'6','₇':'7','₈':'8','₉':'9',
    'ₜ':'t','ₓ':'x','ₙ':'n','ᵢ':'i','ⱼ':'j',
}

COMBINING = {'\u0304': 'bar', '\u0303': 'tilde', '\u0302': 'hat', '\u0307': 'dot'}

ALL_MATH = set(GREEK) | set(OPERATORS) | set(SUPER) | set(SUB) | set(COMBINING) | {'∞', '√', '·', '−'}


def has_math(text):
    return any(c in ALL_MATH for c in text)


def split_protected(line):
    """Split line into ('text',...), ('code',...), ('latex',...), ('url',...)."""
    segs = []
    i = 0
    buf = []
    
    while i < len(line):
        # Markdown link URLs: ](...)
        if line[i:i+2] == '](' and i > 0:
            buf.append(']')
            if buf:
                # flush buf minus the ] we just added... actually let's handle this:
                pass
            # Find matching closing )
            j = i + 2
            depth = 1
            while j < len(line) and depth > 0:
                if line[j] == '(': depth += 1
                elif line[j] == ')': depth -= 1
                j += 1
            # buf already has ']', add the rest as protected
            segs.append(('text', ''.join(buf)))
            buf = []
            segs.append(('url', '(' + line[i+2:j]))
            i = j
            continue
        
        if line[i] == '`':
            if buf:
                segs.append(('text', ''.join(buf)))
                buf = []
            j = i + 1
            while j < len(line) and line[j] == '`':
                j += 1
            ticks = j - i
            close = line.find('`' * ticks, j)
            if close >= 0:
                segs.append(('code', line[i:close + ticks]))
                i = close + ticks
            else:
                buf.append(line[i])
                i += 1
            continue
        
        if line[i] == '$':
            if buf:
                segs.append(('text', ''.join(buf)))
                buf = []
            if i + 1 < len(line) and line[i+1] == '$':
                close = line.find('$$', i + 2)
                if close >= 0:
                    segs.append(('latex', line[i:close+2]))
                    i = close + 2
                else:
                    buf.append('$')
                    i += 1
            else:
                close = line.find('$', i + 1)
                if close >= 0:
                    segs.append(('latex', line[i:close+1]))
                    i = close + 1
                else:
                    buf.append('$')
                    i += 1
            continue
        
        buf.append(line[i])
        i += 1
    
    if buf:
        segs.append(('text', ''.join(buf)))
    
    return segs


def convert_inline(text):
    """Convert Unicode in inline text, each symbol wrapped in $...$."""
    result = []
    i = 0
    while i < len(text):
        c = text[i]
        
        # Combining diacriticals (char + combining mark)
        if i + 1 < len(text) and text[i+1] in COMBINING:
            cmd = COMBINING[text[i+1]]
            result.append(f'$\\{cmd}{{{c}}}$')
            i += 2
            continue
        
        if c in GREEK:
            # Merge Greek + following single letter as one math expression (Δx → $\Delta x$)
            if i + 1 < len(text) and text[i+1].isalpha() and (i + 2 >= len(text) or not text[i+2].isalpha()):
                result.append(f'${GREEK[c]} {text[i+1]}$')
                i += 2
            else:
                result.append(f'${GREEK[c]}$')
                i += 1
            continue
        
        if c in OPERATORS:
            result.append(f'${OPERATORS[c]}$')
            i += 1
            continue
        
        # Middle dot: skip in inline (too many false positives with text separators)
        # Will be handled manually in equations
        if c == '·':
            result.append(c)
            i += 1
            continue
        
        if c == '−':
            result.append('$-$')
            i += 1
            continue
        
        if c == '∞':
            if result:
                prev = result[-1]
                if len(prev) == 1 and prev.isalpha():
                    result.pop()
                    result.append(f'${prev}_\\infty$')
                elif prev.startswith('$') and prev.endswith('$'):
                    result.pop()
                    inner = prev[1:-1]
                    result.append(f'${inner}_\\infty$')
                else:
                    result.append(r'$\infty$')
            else:
                result.append(r'$\infty$')
            i += 1
            continue
        
        if c == '√':
            if i + 1 < len(text) and text[i+1] == '(':
                depth = 0
                j = i + 1
                while j < len(text):
                    if text[j] == '(': depth += 1
                    elif text[j] == ')':
                        depth -= 1
                        if depth == 0: break
                    j += 1
                inner = text[i+2:j]
                result.append(f'$\\sqrt{{{inner}}}$')
                i = j + 1
            else:
                result.append(r'$\sqrt{}$')
                i += 1
            continue
        
        if c in SUPER:
            digits = []
            while i < len(text) and text[i] in SUPER:
                digits.append(SUPER[text[i]])
                i += 1
            s = ''.join(digits)
            result.append(f'$^{{{s}}}$' if len(s) > 1 else f'$^{s}$')
            continue
        
        if c in SUB:
            digits = []
            while i < len(text) and text[i] in SUB:
                digits.append(SUB[text[i]])
                i += 1
            s = ''.join(digits)
            result.append(f'$_{{{s}}}$' if len(s) > 1 else f'$_{s}$')
            continue
        
        result.append(c)
        i += 1
    
    return ''.join(result)


def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    lines = content.split('\n')
    result = []
    in_code_block = False
    
    for line in lines:
        stripped = line.lstrip()
        
        if stripped.startswith('```'):
            in_code_block = not in_code_block
            result.append(line)
            continue
        
        if in_code_block:
            result.append(line)
            continue
        
        if not has_math(line):
            result.append(line)
            continue
        
        # Process inline text only
        segs = split_protected(line)
        converted_parts = []
        for seg_type, seg_text in segs:
            if seg_type in ('code', 'latex', 'url'):
                converted_parts.append(seg_text)
            else:
                converted_parts.append(convert_inline(seg_text))
        
        result.append(''.join(converted_parts))
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write('\n'.join(result))
    
    print(f"Processed: {filepath}")


if __name__ == '__main__':
    for fp in sys.argv[1:]:
        process_file(fp)
