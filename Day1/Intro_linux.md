# Intro to Linux/bash 

*by Maria Eleonora Rossi*

During the course we will use the command line for our exercises. The command line is a text-based interface used widely in bioinformatics to interact with computers.
A bit of familiarity with the command line will help you greatly in our practical; however, here is a little refresher on the commands that are going to help you navigate folders, read files, and launch analyses.

We start the tutorial in the Terminal. Open the Terminal on your machine. You'll likely see something like this
```sh
(base) eleonora@mac ~ % 
```
The `%` is called the prompt; on other machines (e.g. Windows) can also look like this `>`.

## Running commands

Once opened the terminal we are ready to run commands. The command is the first you type, normally can be followed by flags or settings introduced by the symbol `-`, that change the default behaviour of the command.

For example:
```sh
(base) eleonora@mac ~ % ls -l 
total 2536
drwxr-xr-x   4 eleonora  staff      128 26 May  2025 Aggresso
drwx------   2 eleonora  staff       64  8 Jul  2024 Apps
```
The first line after the prompt is the command we are executing. In this case `ls` is listing all files and directories in the current working directory. The flag `-l` forces the `ls` command to print the long version of the output with extended information about the directory's content.

We can also specify to list a specific directory to inspect by adding the name of what we want to inspect after the `-l`
```sh
(base) eleonora@mac ~ %  ls -l Conferences/CGUE_eukaryotics 
total 0
drwxr-xr-x  10 eleonora  staff    320 12 Sep  2024 Expenses claim
-rw-rw----@  1 eleonora  staff  57498 20 Aug  2024 Request-for-approval-to-Travel-Oct-2023_CGUE.docx
```

Whenever we want to know how a command works we can type `man` followed by the command. 
```sh
(base) eleonora@mac ~ %  man ls 
```
Another way is to type `-h` or `--help` after the command or program we want to execute.


## Navigating and working in the filesystem

The nested hierarchy of folders and files on your computer is called the filesystem. When we use the computer, we are used to click on the folder or the file we want to open. On the terminal, we cannot do that with a click, but we need to navigate and access the folder (we normally refer to them as directories) we want to open with specific commands. 

All directories are nested within a previous directory, and you can also have multiple directories in the same location, as you would have, for example, different folders in your Finder app. 

To know where we are in the filesystem of our computer, we can type in the directory `pwd`, short for print working directory.

```sh
(base) eleonora@mac Desktop % pwd
/Users/eleonora/Desktop
```
`/Users/eleonora` is the pathway, and `Desktop` is the directory I am in. Each directory is separated by `/` symbol.

Let's say we want to create a new directory, normally we would do right click -> new. Here we use the command `mkdir` (make directory).

```sh
(base) eleonora@mac Desktop % mkdir Linux_tutorial
```
To check if we have created the directory, we type `ls` (list), which allows us to see all the files and directories present in the directory we're in.

```sh
(base) eleonora@mac Desktop % ls
Linux_tutorial
```
Now we want to go into the directory we just created. We can use `cd` (change directory) and the name of the directory

```sh
(base) eleonora@mac Desktop % cd Linux_tutorial 
(base) eleonora@mac Linux_tutorial % 

```
Type again `pwd` and check the pathway you're in. Is it similar to the previous pwd?

If we want to exit the directory we are in, we need to 'move back' in the folder. The wildcard to do is `..`, the full command is `cd ..`. This command always moves you one directory back. If you want to move 2 directories back, you just need to type `cd ../../`. Remember to delimit each directory with the `/`. 

Type `cd ..`
What is the path now?

Let's create another directory in `Linux_tutorial` and call it `data`. What would you type?
Enter the directory and check the absolute path.

Ok now we have the basics to navigate the filesystem. Let's see how to work with files.

## Working with files

Most times, we need to visualise the files. To do so, there are different commands we can use. `less` is very useful; it loads the file you want to visualise, and you can "scroll" the file by hitting the space bar or actually scrolling with the mouse. The advantage of `less` is that it loads in the memory temporarily only the section of the file you are looking at, meaning that it's handy to look at very large files. 


Two other useful commands are `head` and `tail` 
As the names might suggest, `head` shows you the first part of the file. With the flag `-n`, you can specify how many lines you want to see of the file. 
E.g. if you type `head -n 30 <file-name>` you will visualise the first 30 lines of that file.

A similar concept can be applied to `tail`, where it allows you to check the last part of the file. 
What would you type to check the last 20 lines of the `chicken.fa` file?


## Editing files

The are several text editors you can use to edit a file or a script. One of the easiest to use is `nano`. You can start edit the file by just typing the command and the name of the file you want to edit, or if you want to write a new file just type nano followd by the name you want to give to your file.

```
nano <name-file>
```
Once opened, you can move around with the cursor, write, paste, copy, etc, etc. Once done, you can ext `nano` with `ctrl + x`. If the file has been modified, `nano` will ask if you want to save it, and you say yes by typing `y`. The program will then quit, and you will find the file in the directory. 

## Wildcards and Pattern Search

Vwery often in bioinformatics, we need to run the same command hundreds of times on hundred of files. Wildcards are special characters used to match filenames or text patterns in Unix.

Common wildcards:

- `*` → matches any number of characters
- `?` → matches a single character
- `[]` → matches one character from a set

## Examples that you could use in one of your local folder on your machine

List all `.fasta` files:

```bash
ls *.fasta
```

Match files like `data1.csv` or `dataA.csv`:

```bash
ls data?.csv
```

Match `sample1.fa`, `sample2.fa`, and `sample3.fa`:

```bash
ls sample[123].fa
```

Pattern searching inside files is commonly done with `grep`.

Find lines containing the word “gene”:

```bash
grep "gene" annotations.gff
```

Ignore uppercase/lowercase:

```bash
grep -i "atg" sequences.fa
```

Show line numbers:

```bash
grep -n "exon" genes.gff
```

Search for “ATPase” in all `.txt` files:

```bash
grep "ATPase" *.txt
```

---

# Piping and Redirection

Unix allows commands to work together using pipes and redirection.

## Redirection

`>` sends output to a file:

```bash
ls > files.txt
```

`>>` appends output instead of overwriting:

```bash
echo "new line" >> notes.txt
```

`<` uses a file as input:

```bash
sort < names.txt
```

## Pipes

The pipe symbol `|` sends the output of one command into another command.

Count how many files are in the directory:

```bash
ls | wc -l
```

Show the first matching lines containing “gene”:

```bash
grep "gene" annotations.gff | head
```

Count FASTA sequences by counting headers:

```bash
cat proteins.fa | grep ">" | wc -l
```
Or 
```bash
grep -c ">" proteins.fa 



## Useful resources
These are just some very basic commands and instructions on how to use the command line. Please come back to this tutorial if you need help while running the phylogenomics practical.

For a more in-depth read, a very good starting book for your first steps in bioinformatics is *Practical Computing for Biologists by Haddoc and Dunn*. You can find more info [here](https://practicalcomputing.org/index.html). 

Finally, here is the [cheat sheet](https://linuxstans.com/bash-cheat-sheet/) for some basic bash commands.
 
Good luck and have fun!