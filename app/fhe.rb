module X
  class FHE
    attr_accessor :p1, :p2, :p3, :p4, :q, :k, :e

    def self.random_number(bits)
      OpenSSL::BN::rand(bits).to_i
    end

    def self.random_prime(bits)
      OpenSSL::BN::generate_prime(bits).to_i
    end

    def initialize(gamma,lambda)
      p1_bits = gamma * lambda
      p2_bits = lambda / 2
      p3_bits = lambda - 1
      p4_bits = lambda - 1

      @p1 = self.class.random_prime(p1_bits)
      @p2 = self.class.random_prime(p2_bits)

      # @q = @p1 * @p2 * @p3 * @p4

      ##
      valid_key = false
      while valid_key == false
        begin
          @p3 = self.class.random_prime(p3_bits)
          @p4 = self.class.random_prime(p4_bits)

          @q = @p1 * @p2 * @p3 * @p4

          k1 = self.class.hensel_packing(@p3, @p2, @q)
          k2 = self.class.hensel_packing(@p4, @p2, @q)
          @k = k1.gp(k2)

          @k.inverse

          valid_key = true
        rescue => e
          # puts "*"
        end
      end
      ##

      @e = Multivector2Dff.new [
        self.class.hensel_encoding(Rational(1,2), @p2),
        self.class.hensel_encoding(Rational(1,2), @p2),
        self.class.hensel_encoding(Rational(1,2), @p2),
        self.class.hensel_encoding(Rational(-1,2), @p2)
      ], @q
    end

    def self.hensel_encoding(n, p)
      Xp.new([p], n.numerator, n.denominator).to_i
    end

    def self.hensel_decoding(n, p)
      Xp.new([p], n.numerator, n.denominator).to_r
    end

    def self.hensel_packing(m, p, q)
      p_bits = p.bit_length

      rs = Array.new(4){ FHE.random_number(p_bits / 8) }

      d1 = Rational(rs[0],rs[1])
      d2 = Rational(rs[2],rs[3])

      data = Array.new(3){
        [0,1].sample == 1 ? 1 * FHE.random_number(p_bits / 8) : (-1) * FHE.random_number(p_bits / 8)
      }

      data << (m - data.inject(:+))

      mm = Multivector2D.new data

      data_ = [
        data[0] + d1,
        data[1] - d1,
        data[2] + d2,
        data[3] - d2
      ]

      mm_ = Multivector2D.new data_

      data_p = data_.map{|d|
        hensel_encoding(d, p)
      }

      mm_p = Multivector2Dff.new data_p, q

      mm_p
    end

    def self.hensel_unpacking(mm, p)
      # data = mm.data.map{|d| d % sk.q}

      data_ = mm.data.map{|d|
        hensel_decoding(d, p)
      }

      data_.inject(:+)
    end

    def self.modular_pow( base, power, mod )
      res = 1
      while power > 0
        res = (res * base) % mod if power & 1 == 1
        base = base ** 2 % mod
        power >>= 1
      end
      res
    end

    def self.multivector_to_n(mm,n)
      Multivector2Dff.new [
        modular_pow(mm.e0,n,mm.modulus),
        modular_pow(mm.e1,n,mm.modulus),
        modular_pow(mm.e2,n,mm.modulus),
        modular_pow(mm.e12,n,mm.modulus)
      ], mm.modulus
    end

    def self.undo_multivector_to_n(mm,n)
      Multivector2Dff.new mm.data.map{|d| d % n}, mm.modulus
    end

    def self.standard_packing(n,q)
      b = n.bit_length - 1
      c1 = random_number(b)
      c2 = random_number(b)
      c3 = random_number(b)
      c4 = ((n - (c1 + c2 + c3))) % q

      Multivector2Dff.new [c1,c2,c3,c4], q
    end

    def gp_encrypt(mm)
      k.gp(mm).gp(k.inverse)
    end

    def gp_decrypt(c)
      k.inverse.gp(c).gp(k)
    end

    def encrypt(m)
      mm = self.class.hensel_packing(m,p2,q)
      me = mm.gp(e)
      c_ = self.class.multivector_to_n(me,p1)
      c__ = c_.gp(e)

      c = gp_encrypt(c__)

      c
    end

    def decrypt(c)
      c__ = gp_decrypt(c)

      me = self.class.undo_multivector_to_n(c__,p1)

      m = self.class.hensel_unpacking(me, p2)
      m
    end

    def add(c1,c2)
      c1.add(c2)
    end

    def mul(c1,c2)
      c1.gp(c2)
    end

  end
end
