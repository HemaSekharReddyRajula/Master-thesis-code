import random
def main():
	print "Guess a number between 1 nad 100"
	#randomNumber = 35
	randomNumber = random.randint(1,100)
	found = False
	
	while not found:
		userGuess = input("Your guess: ")
		if userGuess == randomNumber:
			print "You got it"
			found = True
		elif userGuess > randomNumber:
			print "Guess lower"
		else:
			print "Guess higher"
	
if __name__ == "__main__":
	main()
