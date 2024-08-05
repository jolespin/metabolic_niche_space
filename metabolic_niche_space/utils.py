# -*- coding: utf-8 -*-

# Split rule characters
def rule_splitter(rule:str, split_characters:list) -> set:
    """
    Split rule by characters

    Args:
        rule (str): Boolean logical string
        split_characters (list): List of characters to split in rule

    Returns:
        set: Unique symbols in a rule
    """
    # Split the rule into individual symbols
    rule_decomposed = str(rule)
    if split_characters:
        for character in split_characters:
            character = character.strip()
            rule_decomposed = rule_decomposed.replace(character, "")
    unique_symbols = set(rule_decomposed.split())
    return unique_symbols

# Function to parse and evaluate rules
def evaluate_rule(rule:str, items:set, replace={"+":" & ", ",":" | "}) -> bool: # Replace `+` with `&` (AND) and `,` with `|` (OR)
    """
    Evaluate a string of boolean logicals

    Args:
        rule (str): Boolean logical string
        items (set): List of items in rule
        replace (dict, optional): Replace boolean characters. Defaults to {"+":" & ", ",":" | "}.

    Returns:
        bool: Evaluated rule
        
    Usage: 
    # Define the rules
    rules = [
        'R00351', 
        'R01325+R01900,R01324', 
        'R01899+R00268,R00267,R00709', 
        'R00621+R03316,R01700', 
        'R02570', 
        'R07618', 
        'R01197', 
        'R00405,R00432,R00727,R10343', 
        'R02164', 
        'R01082', 
        'R00342,R00361'
    ]

    # List of items to check against
    items = {
        'R00267','R00342','R00361','R00405','R00432',
        'R00621','R03316','R00709','R00727','R01082',
        'R01197','R01900','R01899','R00268','R02164',
        'R02570','R07618','R10343',
    }
    
    # Evaluate and filter the rules
    rule_to_bool = {rule:evaluate_rule(rule, items) for rule in rules}

    rule_to_bool
    # {'R00351': False,
    #  'R01325+R01900,R01324': False,
    #  'R01899+R00268,R00267,R00709': True,
    #  'R00621+R03316,R01700': True,
    #  'R02570': True,
    #  'R07618': True,
    #  'R01197': True,
    #  'R00405,R00432,R00727,R10343': True,
    #  'R02164': True,
    #  'R01082': True,
    #  'R00342,R00361': True}
    """
    # Replace characters for standard logical formatting
    if replace:
        for character_before, character_after in replace.items():
            rule = rule.replace(character_before, character_after)
    
    # Split the rule into individual symbols
    unique_symbols = rule_splitter(rule, replace.values())
    
    # Create a dictionary with the presence of each symbol in the items
    item_to_bool = {sym: (sym in items) for sym in unique_symbols}
    # Parse the rule
    expr = eval(rule, {"__builtins__": None}, item_to_bool)
    # Evaluate the expression
    return expr
