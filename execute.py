import sys

import tkinter as tk
from tkinter import scrolledtext as st
from tkinter import filedialog as fd
from tkinter import messagebox as mb
from tkinter import Button

import Bio as Bio
from Bio.Seq import Seq
from Bio.SeqUtils import GC


seq = None
seq_complementary = None
seq_specific_gene = None
complete_genome = None
genome_20_pb = None
percent_gc_20 = None
translate_seq_list = None

complete_genome_mutation = None
genome_mutation = None
translate_mutation_seq_list = None
list_codon_mutation = []
seq_mutation_complementary = None
list_codon_mutation_adn = None
seq_mutation_specific_gene = None
list_codon_mutation_specific = None

dict_translate = {
    'UUU': 'Phe',
    'UUC': 'Phe',
    'UUA': 'Leu',
    'UUG': 'Leu',

    'CUU': 'Leu',
    'CUC': 'Leu',
    'CUA': 'Leu',
    'CUG': 'Leu',

    'AUU': 'lle',
    'AUC': 'lle',
    'AUA': 'lle',
    'AUG': 'Met',

    'GUU': 'Val',
    'GUC': 'Val',
    'GUA': 'Val',
    'GUG': 'Val',

    'UCU': 'Ser',
    'UCC': 'Ser',
    'UCA': 'Ser',
    'UCG': 'Ser',

    'CCU': 'Pro',
    'CCC': 'Pro',
    'CCA': 'Pro',
    'CCG': 'Pro',

    'ACU': 'Thr',
    'ACC': 'Thr',
    'ACA': 'Thr',
    'ACG': 'Thr',

    'GCU': 'Ala',
    'GCC': 'Ala',
    'GCA': 'Ala',
    'GCG': 'Ala',

    'UAU': 'Tyr',
    'UAC': 'Tyr',
    # 'UAA': 'Stop',
    # 'UAG': 'Stop',

    'CAU': 'His',
    'CAC': 'His',
    'CAA': 'Gin',
    'CAG': 'Gin',

    'AAU': 'Asn',
    'AAC': 'Asn',
    'AAA': 'Lys',
    'AAG': 'Lys',

    'GAU': 'Asp',
    'GAC': 'Asp',
    'GAA': 'Glu',
    'GAG': 'Glu',

    'UGU': 'Cys',
    'UGC': 'Cys',
    # 'UGA': 'Stop',
    'UGG': 'Trp',

    'CGU': 'Arg',
    'CGC': 'Arg',
    'CGA': 'Arg',
    'CGG': 'Arg',

    'AGU': 'Ser',
    'AGC': 'Ser',
    'AGA': 'Arg',
    'AGG': 'Arg',

    'GGU': 'Gly',
    'GGC': 'Gly',
    'GGA': 'Gly',
    'GGG': 'Gly',
}

list_codon_origin = None
list_codon_adn = None
list_codon_specific = None

class App:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title('Algoritmo Dogma Central')
        self.add_menu()
        self.scrolledtext1 = st.ScrolledText(self.root, width=120, height=30)
        self.scrolledtext1.grid(column=0,row=0, padx=10, pady=10)        
        self.save_button = Button(self.root, text="Guardar mutación", width=15, command=self.generate_mutation)
        
        self.root.mainloop()


    def add_menu(self):
        menubar1 = tk.Menu(self.root)
        self.root.config(menu=menubar1)
        choice1 = tk.Menu(menubar1, tearoff=0)
        choice1.add_command(label='Guardar archivo', command=self.save)
        choice1.add_command(label='importar archivo', command=self.import_file)
        choice1.add_separator()
        choice1.add_command(label='Visualizar Genoma Completo', command=self.get_complete_genome)
        choice1.add_command(label='ADN Complementario', command=self.adn_complementary)
        choice1.add_command(label='Transcripción (gen específico)', command=self.specific_gene)
        choice1.add_command(label='Contenido de G/C', command=self.percent_gc)
        choice1.add_command(label='Traducción', command=self.translate)
        choice1.add_command(label='Generar mutaciones en un gen específico', command= self.mutation)
        choice1.add_command(label='Identificación de la cadena peptídica generada', command=self.identification_pep)
        choice1.add_separator()
        choice1.add_command(label='Salir', command=self.exit_app)
        menubar1.add_cascade(label='Archivo', menu=choice1) 


    def show_info(self, text=''):
        self.scrolledtext1.delete('1.0', tk.END) 
        self.scrolledtext1.insert('1.0', text)


    def show_popup(self, title='Error', text='', delete_window=False):
        if delete_window:
            self.scrolledtext1.delete('1.0', tk.END) 
        mb.showinfo(title, text)


    def get_complete_genome(self):
        global list_codon_origin
        if not list_codon_origin:
            self.show_popup(text='Por favor importe su archivo con el genoma completo')
        self.show_info(text=list_codon_origin)


    def adn_complementary(self):
        global seq, seq_complementary, complete_genome, list_codon_adn

        if not complete_genome:
            self.show_popup(text='Por favor importe su archivo con el genoma completo')
            return None

        if not seq_complementary:
            seq_complementary = seq.complement()
            list_codon_adn = [seq_complementary[i:i+3] for i in range(0, len(seq_complementary), 3)]
        self.show_info(text=list_codon_adn[:20])


    def specific_gene(self):
        global seq, seq_specific_gene, list_codon_specific

        if not seq:
            self.show_popup(text='Por favor importe su archivo con el genoma completo')
            return None
        seq_specific_gene = seq.transcribe()
        list_codon_specific = [seq_specific_gene[i:i+3] for i in range(0, len(seq_specific_gene), 3)]
        self.show_info(text=list_codon_specific[:20])


    def percent_gc(self):
        global seq_specific_gene, percent_gc_20

        if not seq_specific_gene:
            self.show_popup(text='Por favor genere su transcripción (gen específico)')
            return None

        percent_gc_20 = GC(seq_specific_gene[:60])
        self.show_info(text=percent_gc_20)


    def translate(self):
        global seq_specific_gene, translate_seq_list, dict_translate, list_codon_specific

        if not seq_specific_gene:
            self.show_popup(text='Por favor genere su transcripción (gen específico)')
            return None

        translate_seq_list = []
        init_position = list_codon_specific.index('AUG')
        for i in list_codon_specific[init_position:]:
            translate = dict_translate.get(str(i)) 
            if translate: 
                translate_seq_list.append([str(i), translate])
                continue
            break
        self.show_info(text=translate_seq_list[:20])


    def mutation(self):
        global list_codon_origin

        if not list_codon_origin:
            self.show_popup(text='Por favor importe su archivo con el genoma completo')
            return None

        show_list_codon_origin = ''
        count = 0
        for codon in list_codon_origin:
            show_list_codon_origin += '{}. {}\n'.format(count, codon)
            count += 1
        self.show_popup(title='Intrucciones', text='''
1. Para editar: Encuentre el número del codon a editar y solo modifique las letras en mayusculas.
2. Para agregar: Solo añada un número con un punto y espacio; solo 3 letras en mayusculas.
3. Para eliminar: Elimine la linea completa (número y letras) que desea.
        ''')
        self.show_info(text=show_list_codon_origin)

        self.save_button.grid(row=4, column=0, columnspan=3, padx=5, pady=5)


    def generate_mutation(self):
        global list_codon_origin, list_codon_mutation, complete_genome_mutation

        if not list_codon_origin:
            self.show_popup(text='Por favor importe su archivo con el genoma completo')
            return None

        text = self.scrolledtext1.get('1.0', tk.END)
        self.scrolledtext1.delete('1.0', tk.END)
        self.save_button.grid_remove()
        self.scrolledtext1.insert('1.0', list_codon_origin[:20])
        self.scrolledtext1.insert(tk.END, '\n \n' )
        self.scrolledtext1.pack()
        self.scrolledtext1.tag_config('modified', foreground='red')
        position = 0

        for element in text.split('\n'):
            if element != '':
                y = element.split('. ')[1].strip().upper()
                list_codon_mutation.append(y)
                if position < 20:
                    try:
                        if y != list_codon_origin[position]:
                            y = self.scrolledtext1.insert(tk.END, y + ' ', 'modified')
                        else:
                            y = self.scrolledtext1.insert(tk.END, y + ' ')
                    except IndexError:
                        y = self.scrolledtext1.insert(tk.END, y + ' ', 'modified')
                position += 1

                
        complete_genome_mutation = Seq(''.join(list_codon_mutation))


    def identification_pep(self):
        global list_codon_origin, translate_seq_list, list_codon_mutation, complete_genome_mutation, dict_translate, seq_mutation_complementary, seq_mutation_specific_gene, list_codon_mutation_specific
        
        if not list_codon_origin:
            self.show_popup(text='Por favor importe su archivo con el genoma completo')
            return None
        elif not translate_seq_list:
            self.show_popup(text='Por favor genere la traducción de la cadena original')
            return None
        elif not list_codon_mutation:
            self.show_popup(text='Por favor genere una mutación')
            return None

        if not seq_mutation_complementary:
            seq_mutation_complementary = complete_genome_mutation.complement()
            list_codon_mutation_adn = [seq_mutation_complementary[i:i+3] for i in range(0, len(seq_mutation_complementary), 3)]
        if not seq_mutation_specific_gene:
            seq_mutation_specific_gene = complete_genome_mutation.transcribe()
            list_codon_mutation_specific = [seq_mutation_specific_gene[i:i+3] for i in range(0, len(seq_mutation_specific_gene), 3)]

        translate_mutation_seq_list = []
        init_position = list_codon_mutation_specific.index('AUG')
        for i in list_codon_mutation_specific[init_position:]:
            translate = dict_translate.get(str(i))
            if translate:
                translate_mutation_seq_list.append([str(i), translate])
                continue
            break
            
        self.scrolledtext1.delete('1.0', tk.END)
        self.scrolledtext1.insert('1.0', translate_seq_list[:20])
        self.scrolledtext1.insert(tk.END, '\n \n' )
        self.scrolledtext1.pack()
        self.scrolledtext1.tag_config('modified', foreground='red')
        position = 0

        for element in translate_mutation_seq_list:
            if element != '':
                if position < 20:
                    try:
                        if element != translate_seq_list[position]:
                            element = self.scrolledtext1.insert(tk.END, str(element) + ' ', 'modified')
                        else:
                            element = self.scrolledtext1.insert(tk.END, str(element) + ' ')
                    except IndexError:
                        element = self.scrolledtext1.insert(tk.END, str(element) + ' ', 'modified')
                position += 1


    def exit_app(self):
        sys.exit()


    def save(self):
        file_name = fd.asksaveasfilename(initialdir = '/',title = 'Guardar como',filetypes = (('txt files','*.txt'),('todos los archivos','*.*')))
        if file_name != '':
            archi1 = open(file_name, 'w', encoding='utf-8')
            archi1.write(self.scrolledtext1.get('1.0', tk.END))
            archi1.close()
            self.show_popup('Información', 'Los datos fueron guardados en el archivo.')


    def import_file(self):
        global complete_genome, seq, list_codon_origin, genome_20_pb

        file_name = fd.askopenfilename(initialdir = '/',title = 'Seleccione archivo',filetypes = (('txt files','*.txt'),('todos los archivos','*.*')))
        if file_name != '':
            archi1 = open(file_name, 'r', encoding='utf-8')
            content = archi1.read()
            complete_genome = content.replace('\n', '').replace(' ', '').split('genome')[1]
            list_codon_origin = [complete_genome[i:i+3] for i in range(0, len(complete_genome), 3)]
            genome_20_pb = ''.join(list_codon_origin[:20])
            archi1.close()
            seq = Seq(complete_genome)
            self.show_info(text=list_codon_origin)


def main():
    my_app = App()
    return 0

if __name__ == '__main__':
    main()